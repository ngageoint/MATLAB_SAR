function [ readerobj ] = open_cphd30_reader( filename )
%OPEN_CPHD30_READER Intiates a reader object for Compensated Phase
%History Data (CPHD) version 3.0.  Metadata is returned in CPHD "X" format
%(newer than version 3.0).
%
% Caveats: 1) Only works on CPHD version 3.0.
%          2) Reader currently doesn't handle multiple channels
%          3) Assumes an order to the pulses, which while reasonable is not
%             guaranteed in CPHD spec.
%          4) The AmpSF0 scaling factor is already applied to WBVECTORS
%             upon reading, so there is no need to "reapply" it outside the
%             calling of this function.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Read metadata and process for later processing
cphd_preamble=read_cphd_preamble(filename);
if isfield(cphd_preamble,'Npulses') % For older versions of CPHD (< 3.0)
    cphd_preamble.Nvectors=cphd_preamble.Npulses; % "Vectors" field naming convention
end
switch cphd_preamble.PHDataType
    case 'cmplxn'
        matlab_datatype='bit4';
        datatype_bytes=0.5;
    case 'cmplxb'
        matlab_datatype='int8';
        datatype_bytes=1;
    case 'cmplxs'
        matlab_datatype='int16';
        datatype_bytes=2;
    case 'cmplxf'
        matlab_datatype='single';
        datatype_bytes=4;
    otherwise
        error('OPEN_CPHD30_READER:UNRECOGNIZED_DATATYPE','Unrecognized data type.');
end
% Adjust frequencies to be true, not offset values, if reference frequency
% is available
if isfield(cphd_preamble,'FreqReferenceIndex')&&cphd_preamble.FreqReferenceIndex&&exist('cphd_ref_freq','file')
    ref_freq=cphd_ref_freq();
    if isfield(cphd_preamble,'NominalCenterFreq')
        cphd_preamble.NominalCenterFreq = cphd_preamble.NominalCenterFreq + ref_freq;
    end
    cphd_preamble.FreqReferenceIndex=0;
else
    ref_freq=0;
end
% Open file to prepare for reading
fid=fopen(filename,'r',cphd_preamble.endian);

% Setup reader object
vb_fieldnames_native = fieldnames(cphd_preamble.VectorParameters);
vb_fieldnames_cphdx = convert_vb_fieldnames_to_cphdx(vb_fieldnames_native);
if cphd_preamble.Interleaved % More generic read function
    [unused, all_nbdata] = read_data(1:cphd_preamble.Nvectors, [], 1);
else % This way is much faster though.
    all_nbdata = read_all_noninterleaved_nb();
end
all_nbdata.SC0 = all_nbdata.SC0 + ref_freq; % Adjust for reference frequency
meta=meta2cphdx_cphd30(cphd_preamble, all_nbdata); % Convert metadata to CPHDX format
meta.native.cphd30 = cphd_preamble; % Save original format

readerobj.read_cphd=@read_data;
readerobj.get_meta=@() meta;
readerobj.close=@() fclose(fid);

    function [wbvectors, nbdata] = read_data(pulse_indices, sample_indices, channels)
        % Parse input parameters
        if (nargin<1)||strcmpi(pulse_indices,'all')
            pulse_indices=1:cphd_preamble.Nvectors;
        end
        if (nargin<2)||strcmpi(sample_indices,'all')
            sample_indices=1:cphd_preamble.Nsamples;
        end
        if (nargin<3)||strcmpi(channels,'all')
            channels=1:cphd_preamble.Nchannels;
        end

        % Setup variables
        if cphd_preamble.Interleaved
            % Preallocate nbdata datastructure
            for i = 1:length(vb_fieldnames_cphdx)
                num_vb_elements = ceil(cphd_preamble.VectorParameters.(vb_fieldnames_native{i})/8);
                nbdata.(vb_fieldnames_cphdx{i}) = zeros(length(pulse_indices),num_vb_elements,length(channels));
            end
        end
        [nb_offsets, wb_offsets]=pulse_offsets(pulse_indices, channels); % Offsets in file to each pulse
        wbvectors=zeros(length(sample_indices),length(pulse_indices),length(channels),'single'); % For storing binary pulse data
        for i=1:length(channels)
            if cphd_preamble.Interleaved
                for j=1:length(pulse_indices)
                    fseek(fid, nb_offsets(j,i),'bof');
                    read_nb();
                    read_wb(); % In interleaved mode, wideband starts immediately after narrowband.  No need to fseek.
                end
            else
                if nargout>1
                    for fn_index = 1:numel(vb_fieldnames_cphdx)
                        nbdata.(vb_fieldnames_cphdx{fn_index}) =...
                            all_nbdata.(vb_fieldnames_cphdx{fn_index})(pulse_indices,:,i);
                    end
                end
                for j=1:length(pulse_indices)
                    fseek(fid, wb_offsets(j,i),'bof');
                    read_wb();
                end
            end
        end
        
        % Nested functions (read_nb, read_wb, pulse_offsets)
        function read_nb
            nbdata.ChannelNumber(j,i)  = fread(fid, 1, '*uint32')+1; % Zero-based, so add one
            nbdata.VectorNumber(j,i)   = fread(fid, 1, '*uint32')+1; % Zero-based, so add one
            for k = 3:length(vb_fieldnames_cphdx)
                nbdata.(vb_fieldnames_cphdx{k})(j,:,i) =...
                    fread(fid, cphd_preamble.VectorParameters.(vb_fieldnames_native{k})/8, '*double');
            end
        end

        function read_wb
            if ~isempty(sample_indices)
                ChannelNumber        = fread(fid, 1, '*uint32')+1; % Zero-based, so add one
                VectorNumber         = fread(fid, 1, '*uint32')+1; % Zero-based, so add one
                temp                 = fread(fid, 2*cphd_preamble.Nsamples, ['*' matlab_datatype]);
                temp                 = complex(temp(1:2:2*cphd_preamble.Nsamples), temp(2:2:2*cphd_preamble.Nsamples));
                temp                 = single(temp(sample_indices));
                if any(strncmpi(vb_fieldnames_cphdx,'AmpSF',5))
                    wbvectors(:,j,i) = all_nbdata.AmpSF(pulse_indices(j),i) * temp; % Apply scaling factor
                else
                    wbvectors(:,j,i) = temp;
                end
            end
        end

        % Computes offset for a given pulse from the beginning of the the CPHD
        % file (cphd_preamble.pulseStart field from the read_cphd_preamble function.)
        % Assumes data is written in channel order and pulse number order.
        % Also assumes that all pulses from a single channel are written
        % together (pulse number is more quickly increasing than channel number).
        % ACCORDING TO CPHD SPEC, THE IS NOT GUARANTEED!  THIS WILL NOT WORK IF
        % PULSES ARE NOT STORED IN CONSECUTIVE ORDER!!!
        function [nb_offsets, wb_offsets] = pulse_offsets(pulse_numbers, channel_numbers)
            NB_RECORD_LENGTH=(2*4)+(13*8); % Two uint32 fields, three 3-element doubles, and four single-element doubles
            if any(strncmpi(vb_fieldnames_cphdx,'AmpSF',5))
                NB_RECORD_LENGTH = NB_RECORD_LENGTH + 8;
            end
            WB_RECORD_LENGTH=(2*4)+(2*cphd_preamble.Nsamples*datatype_bytes); % Two uint32 field and Nsamples of complex data
            if cphd_preamble.Interleaved
                nb_offsets=cphd_preamble.pulseStart+((pulse_numbers(:)-1)*(NB_RECORD_LENGTH+WB_RECORD_LENGTH));
                wb_offsets=nb_offsets+NB_RECORD_LENGTH;
            else
                nb_offsets=cphd_preamble.pulseStart+((pulse_numbers(:)-1)*NB_RECORD_LENGTH);
                wb_offsets=cphd_preamble.pulseStart+(cphd_preamble.Nvectors*NB_RECORD_LENGTH)+...
                    ((pulse_numbers(:)-1)*WB_RECORD_LENGTH);
            end
            % TODO: Compute offsets due to channel numbers
        end
    end

    % Assumes a single channel dataset
    function nbdata = read_all_noninterleaved_nb()
        fseek(fid, cphd_preamble.pulseStart, 'bof');
        vb_elements_double = sum(cell2mat(struct2cell(cphd_preamble.VectorParameters)))/8;
        % Read all metadata in as double, since a single fread is much faster.
        nbdata_array = fread(fid, [vb_elements_double cphd_preamble.Nvectors],'*double');
        % Cast channel/vector number fields into correct type
        channel_vector = typecast(nbdata_array(1,:),'uint32').'+1;
        % Convert arrays into structure
        nbdata.ChannelNumber = channel_vector(1:2:end);
        nbdata.VectorNumber = channel_vector(2:2:end);
        index = 2;
        for i = 3:length(vb_fieldnames_cphdx)
            num_vb_elements = cphd_preamble.VectorParameters.(vb_fieldnames_native{i})/8;
            nbdata.(vb_fieldnames_cphdx{i}) = nbdata_array(index+(0:(num_vb_elements-1)),:).';
            index = index + num_vb_elements;
        end
    end

    % Convert vector-based fieldnames from a variety of CPHD versions to a standard set of fieldnames
    function cphdx_fieldnames = convert_vb_fieldnames_to_cphdx(cphd30_fieldnames)
        cphdx_fieldnames = cphd30_fieldnames;
        % Converstion for CPHD 1.2 to CPHD 3.0
        replace_fieldname('Fx1', 'SC0');
        replace_fieldname('RxTime', 'RcvTime');
        replace_fieldname('RxPos', 'RcvPos');
        replace_fieldname('PulseNumber','VectorNumber');
        % Conversion from CPHD 3.0 to CPHDX
        replace_fieldname('SRP', 'SRPPos');
        replace_fieldname('AmpSF0', 'AmpSF');
        replace_fieldname('FxStepSize', 'SCSS');
        
        function replace_fieldname(old_fieldname, new_fieldname)
            indices_found = strcmp(old_fieldname,cphdx_fieldnames);
            if any(indices_found)
                cphdx_fieldnames{indices_found} = new_fieldname;
            end
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////