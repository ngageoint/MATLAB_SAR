function [ readerobj ] = open_crsd_reader( filename )
%OPEN_CRSD_READER Intiates a reader object for Compensated Receive
%Signal Data (CRSD) extended format.
%
% A lot of overlap between this and open_cphd_reader.  Perhaps a cleaner
% way to reuse code between them.
%
% Written by: Wade Schwartzkopf, NGA/Resarch
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% File header
fid = fopen(filename,'r','b','UTF-8'); % All CRSD is big-endian
file_header=read_crsd_file_header(fid);

%% XML metadata
fseek(fid,file_header.XML_BLOCK_BYTE_OFFSET,'bof'); % Should be here already, but just to be sure
xml_string = fread(fid,file_header.XML_BLOCK_SIZE,'*char')';
fseek(fid,2,'cof'); % Pass over \f\n section terminator
xml_meta = sicdxml2struct(xmlread(java.io.StringBufferInputStream(xml_string)));

% This might make it work in some old CPHD tools:
% xml_meta.Data.ArraySize = xml_meta.Data.Channel;
% xml_meta.Data.NumCPHDChannels = xml_meta.Data.NumCRSDChannels;
% xml_meta.GeoData.SCP = xml_meta.ReferenceGeometry.CRP;
% xml_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate = xml_meta.Channel.Parameters.Fs;
% xml_meta.RadarCollection.Waveform.WFParameters.RcvFreqStart = xml_meta.Channel.Parameters.F0Ref;
% xml_meta.CollectionInfo = xml_meta.CollectionID;

%% Read in vector-based metadata
% Currently, we read ALL vector-based metadata upon opening the reader.
% This normally is quick, even for large datasets.  Perhaps at some point
% this should be moved into the read_data function where it would only read
% from file the vector-based metadata for the selected vectors, but it
% hardly seems necessary as fast as this goes.
fseek(fid,file_header.PVP_BLOCK_BYTE_OFFSET,'bof'); % Should be here already, but just to be sure
vectorParametersFlatten = xml_meta.PVP;
if isfield(vectorParametersFlatten,'RcvAntenna')
    vectorParametersFlatten = rmfield(vectorParametersFlatten, 'RcvAntenna');
    vectorParametersFlatten = setstructfields(vectorParametersFlatten, xml_meta.PVP.RcvAntenna);
end
if isfield(vectorParametersFlatten,'TxPulse')
    vectorParametersFlatten = rmfield(vectorParametersFlatten, 'TxPulse');
    vectorParametersFlatten = setstructfields(vectorParametersFlatten, xml_meta.PVP.TxPulse);
end
if isfield(vectorParametersFlatten,'TxAntenna')
    vectorParametersFlatten = rmfield(vectorParametersFlatten, 'TxAntenna');
    vectorParametersFlatten = setstructfields(vectorParametersFlatten, xml_meta.PVP.TxPulse.TxAntenna);
end
vectorParametersCell = fieldnames(vectorParametersFlatten);
% Iterate through channels to extract vector-based metadata
bytes_per_element = 8; % Everything in vector-based metadata is 8-byte type
for i = 1:xml_meta.Data.NumCRSDChannels
    % Read all the vector-based metadata for a single channel.  Each column
    % represents metadata for a unique vector.  Each row represents a
    % single parameter for all vectors.
    vb_array = fread(fid, [(xml_meta.Data.NumBytesPVP/bytes_per_element) ...
        xml_meta.Data.Channel(i).NumVectors], 'double');
    % Distribute vector-based metadata into a structure
    for j = 1:numel(vectorParametersCell)
        vbp_all(i).(vectorParametersCell{j}) = ...
            vb_array(vectorParametersFlatten.(vectorParametersCell{j}).Offset+...
            (1:(vectorParametersFlatten.(vectorParametersCell{j}).Size)),:).';
    end
end

%% Setup readers for wideband vector data
% Convert CRSD type string to MATLAB type
switch xml_meta.Data.SignalArrayFormat
    case 'CI2'
        matlab_datatype='int8';
        datatype_bytes=1;
    case 'CI4'
        matlab_datatype='int16';
        datatype_bytes=2;
    case 'CF8'
        matlab_datatype='single';
        datatype_bytes=4;
    otherwise
        error('OPEN_CRSD_READER:UNRECOGNIZED_DATATYPE','Unrecognized data type.');
end
% CRSD does not guarantee that all channels have the same number of
% vectors and samples, so we can't treat as a simple 3D array of size:
% vectors x samples x channels.
% Calculate starting position of CRSD data for each channel
channel_offsets = zeros(xml_meta.Data.NumCRSDChannels,1) + file_header.SIGNAL_BLOCK_BYTE_OFFSET;
for i = 2:xml_meta.Data.NumCRSDChannels
    channel_offsets(i) = channel_offsets(i-1) + (2*datatype_bytes* ...
        xml_meta.Data.Channel(i-1).NumSamples*xml_meta.Data.Channel(i-1).NumVectors);
end
mm_object = cell(xml_meta.Data.NumCRSDChannels,1); % Preallocate array in case memory mapped IO can be used
% Initialize a reader for each channel
for i = 1:xml_meta.Data.NumCRSDChannels
    try % Try memory-mapped IO first
        mm_object{i} = memmapfile(filename,'Offset',channel_offsets(i),...
            'Format',{matlab_datatype [2 xml_meta.Data.Channel(i).NumSamples xml_meta.Data.Channel(i).NumVectors] 'wbvectors'},...
            'Repeat',1);
        % The following line will give out-of-error memory if memmapfile
        % not possible for this data, and thus resort to using fread.
        mm_object{i}.data.wbvectors(1,1);
        % The following only needs to be done once
        if i==xml_meta.Data.NumCRSDChannels
            read_chip_fun = @chip_with_mm;
            [comptype,varsize,compend]=computer;
            need_to_swap_bytes='B'~=upper(compend);
            readerobj.close=@() clear('mm_object'); % Frees file
            fclose(fid); % Don't need to leave the file open for fread anymore
        end
    catch % If memory-mapping fails (on 32-bit machines with large datasets)
        read_chip_fun = @chip_with_fread;
        fread_type = matlab_datatype;
        if strcmp(matlab_datatype,'single')
            fread_type = 'float32';
        end
        readerobj.close=@() fclose(fid);
    end
end

%% Specify reader object methods-- close() method already defined above
readerobj.read_raw = @read_data;
readerobj.get_meta = @() xml_meta;

    %% READ_RAW method of this reader object
    function [wbvectors, nbdata] = read_data(pulse_indices, sample_indices, channel)
        % Compute default input parameters
        if (nargin<3)
            channel=1;
        end
        if (nargin<1)||strcmpi(pulse_indices,'all')
            pulse_indices=1:double(xml_meta.Data.Channel(channel).NumVectors);
        end
        if (nargin<2)||strcmpi(sample_indices,'all')
            sample_indices=1:double(xml_meta.Data.Channel(channel).NumSamples);
        end

        % Read CRSD data
        if ~isempty(sample_indices)
            wbvectors = read_chip_fun(pulse_indices, sample_indices, channel);
            % Rearrange to be complex, instead of IQ interleaved
            wbvectors = permute(complex(wbvectors(1,:,:), wbvectors(2,:,:)),[2 3 1]);
            % Scale amplitude by AmpSF factor
            if isfield(vbp_all(channel), 'AmpSF')
                wbvectors = bsxfun(@times, vbp_all(channel).AmpSF(pulse_indices).', single(wbvectors));
                % Older versions of MATLAB might not have BSXFUN
                %  for k = 1:numel(pulse_indices)
                %      wbvectors(:,k) = wbvectors(:,k) .* vbp_all(channel).AmpSF(pulse_indices(k));
                %  end
            end
        else
            wbvectors = [];
        end

        % Copy selected vector-based data to a structure
        if nargout>1
            fldnames = fieldnames(vbp_all);
            for fn_index = 1:numel(fldnames)
                nbdata.(fldnames{fn_index}) = vbp_all(channel).(fldnames{fn_index})(pulse_indices,:);
            end
        end
    end

    %% Function for reading data with memory mapped files
    % Faster and easier
    function data_out = chip_with_mm(pulse_indices, sample_indices, channel)
        data_out = mm_object{channel}.data.wbvectors(:,sample_indices,pulse_indices);
        if need_to_swap_bytes, data_out = swapbytes(data_out); end;
    end

    %% Function for reading data with fread
    % Works on even low-memory systems
    function data_out = chip_with_fread(pulse_indices, sample_indices, channel)
        % If consecutive pulses and samples, we can avoid looping
        if all(diff(pulse_indices)==1)&&all(diff(sample_indices)==1)
            fseek_to_ps(pulse_indices(1), sample_indices(1), channel);
            data_out=fread(fid,double([2*length(sample_indices),length(pulse_indices)]),...
                [num2str(2*length(sample_indices)) '*' fread_type '=>' matlab_datatype],...
                2*double((xml_meta.Data.Channel(channel).NumSamples-length(sample_indices))*datatype_bytes));
            data_out = reshape(data_out,[2,length(sample_indices),length(pulse_indices)]);
        else % Arbitrary sets of pulses
            % Loop through all pulses
            data_out = zeros(2, length(sample_indices),length(pulse_indices));
            for k = 1:length(pulse_indices)
                fseek_to_ps(pulse_indices(k), sample_indices(1), channel);
                temp = fread(fid, [2 double(xml_meta.Data.Channel(channel).NumSamples)], ...
                    [fread_type '=>' matlab_datatype]);
                data_out(:,:,k) = temp(:,sample_indices);
            end
        end
        
        % Function to seek to a given pulse/sample
        function fseek_to_ps(pulse, sample, channel)
            if fseek(fid,channel_offsets(channel)+... % Beginning of data
                (double(pulse-1)*(double(xml_meta.Data.Channel(channel).NumSamples)*2*datatype_bytes))+... % Skip to first row of interest
                (double(sample-1)*2*datatype_bytes),'bof') % Skip to first column of interest
                error('OPEN_CRSD_READER:EOF',['Attempt to read past end of file.  '...
                 'File possibly corrupt.']);
            end
        end
    end

end

%% Parse the CRSD file header.
% Just use labels in file header as the structure field names.
function meta = read_crsd_file_header(fid)

fseek(fid,0,'bof');
[crsd_string, version] = strtok(fgets(fid),'/'); % First line is just CRSD/version
meta.VERSION = deblank(version(2:end));
if ~strcmp(crsd_string,'CRSD')||~strncmp(meta.VERSION,'1.0.0',5)
    error('OPEN_CRSD_READER:INVALID_FORMAT','Not a valid CRSD file.');
end
kvps = fgetl(fid); % Key-Value pairs
while ~isequal(kvps,12) % \f\n is section terminator is CRSD
    [key, value] = strtok(kvps, ':=');
    key = strtrim(key); % Remove spaces to make sure its a valid fieldname
    try
        meta.(key) = strtrim(value(3:end));
        doubleval = str2double(meta.(key));
        if(~isnan(doubleval)) % If value is numeric, store as such
            meta.(strtrim(key)) = doubleval;
        end
    catch
        error('OPEN_CRSD_READER:INVALID_FORMAT','Not a valid CRSD file.');
    end
    kvps = fgetl(fid); % Key-Value pairs
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////