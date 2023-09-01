function convert_to_cphdx( input_ph_filename, output_cphd_file, varargin )
%CONVERT_TO_CPHD Converts a phase history dataset to CPHD file format
%
% CONVERT_TO_CPHD(INPUT_PH_FILENAME, OUTPUT_CPHD_FILENAME, 'PropertyName', PropertyValue, ...)
%
%       Property name     Description
%       pulse_indices     Pulses from the input phase history dataset to be
%                         written to the output CPHD file.  Default is all
%                         pulses.
%       channels          Channels from the input phase history dataset to
%                         be written to the output CPHD file.  Default is
%                         all channels.
%       SampleType        Datatype to be written to the output CPHD file.
%                         Default is whatever the sample type was closest
%                         to in the input file.  Valid options are 'CI2',
%                         'CI4', and 'CF8'.
%       crsd              Write raw, rather than compensated, phase history
%                         to a CPHD file.  Of course, this option is only
%                         possible for file formats that actually contain
%                         raw, rather than compensated, phase history.
%                         Default is false.
%
% Potential properties (not currently implemented):
% 1) Could allow for translation between FX and TOA domain types.
% 2) Could also add a feature to select samples, but then one would have to
% recalculate some narrowband data (at least start frequency), and assure
% that selected samples were evenly spaced (or maybe more simply just
% consecutive) since CPHD only supports a single Fx_SS for each pulse.
%
% Author: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
% Open phase history file
ph_reader = open_ph_reader(input_ph_filename);
cphd_meta = ph_reader.get_meta();
p1 = inputParser;
p1.KeepUnmatched=false;
all_vectors = 1:cphd_meta.Data.Channel(1).NumVectors;
p1.addParamValue('pulse_indices', all_vectors, @(x) all(ismember(x,all_vectors)));
all_channels = 1:numel(cphd_meta.Data.Channel);
p1.addParamValue('channels', all_channels, @(x) all(ismember(x,all_channels)));
p1.addParamValue('SampleType', cphd_meta.Data.SignalArrayFormat, @ischar);
% Switching domain type not currently supported, but we leave a placeholder for this option
% p1.addParamValue('DomainType', cphd_meta.Global.DomainType, ...
%     @(x) sum(strcmp(x, {'FX','TOA'})==1));
% Option to writes raw data instead of compensated into CPHD format
p1.addParamValue('crsd', false, @(x) isscalar(x)&&islogical(x)&&isfield(ph_reader,'read_raw'));
p1.FunctionName = mfilename;
p1.parse(varargin{:});
% If we are changing any data, update metadata
if ~any(strcmp('pulse_indices',p1.UsingDefaults))
    for i = all_channels
        cphd_meta.Data.Channel(i).NumVectors = length(p1.Results.pulse_indices);
    end
end
if ~any(strcmp('SampleType',p1.UsingDefaults))
    cphd_meta.Data.SignalArrayFormat = p1.Results.SampleType;
end    
% Convert CPHD type string to MATLAB type string
switch cphd_meta.Data.SignalArrayFormat
    case 'CI2'
        matlab_type = 'int8';
    case 'CI4'
        matlab_type = 'int16';
    case 'CF8'
        matlab_type = 'float32';
    otherwise
        error('CONVERT_TO_CPHDX:UNRECOGNIZED_DATATYPE','Unrecognized data type.');
end
% Are we reading raw or compensated data?
if p1.Results.crsd
    read_function = ph_reader.read_raw;
    format = 'CRSD';
else
    read_function = ph_reader.read_cphd;
    format = 'CPHD';
end

%% Update CPHD XML metadata
% Requested channels
if strcmp(format,'CPHD')
    cphd_meta.Data.NumCPHDChannels = numel(p1.Results.channels);
elseif strcmp(format,'CRSD')
    cphd_meta.Data.NumCRSDChannels = numel(p1.Results.channels);
end
cphd_meta.Data.Channel = cphd_meta.Data.Channel(p1.Results.channels);
cphd_meta.Channel = cphd_meta.Channel(p1.Results.channels);
% Populate PVP data if left off
if ~isfield(cphd_meta,'PVP')
    [wb,pvp]=read_function(p1.Results.pulse_indices(1),'all',p1.Results.channels(1));
    fn = fieldnames(pvp);
    pvp_offset = 0;
    for i = 1:numel(fn)
        pvp_size = numel(pvp.(fn{i}));
        if pvp_size == 1
            if strcmp(fn{i}, 'SIGNAL')
                pvp_format = 'I8';
            else
                pvp_format = 'F8';
            end
        elseif pvp_size == 3
            if any(strcmp(fn{i}, {'RcvEB','TxEB'}))
                pvp_format = 'DCX=F8;DCY=F8';
            elseif strcmp(fn{i}, '')
                pvp_format = 'PhiXC=F8;FxC=F8;FxRate=F8';
            else
                pvp_format = 'X=F8;Y=F8;Z=F8';
            end
        end
        pvp_struct = struct('Offset',pvp_offset,'Size',pvp_size,'Format',pvp_format);
        if strcmp(format, 'CRSD')
            switch fn{i}
                case {'RcvACX','RcvACY','RcvEB'}
                    cphd_meta.PVP.RcvAntenna.(fn{i})= pvp_struct;
                case {'TxTime','TxPos','TxVel','FX1','FX2','TXmt','TxLFM'}
                    cphd_meta.PVP.TxPulse.(fn{i})= pvp_struct;
                case {'TxACX','TxACY','TxEB'}
                    cphd_meta.PVP.TxPulse.TxAntenna.(fn{i})= pvp_struct;
                otherwise
                    cphd_meta.PVP.(fn{i})= pvp_struct;
            end
        else
            cphd_meta.PVP.(fn{i})= pvp_struct;
        end
        pvp_offset = pvp_offset + pvp_size;
    end
    cphd_meta.Data.NumBytesPVP = pvp_offset * 8;
end
% Integer types need appropriate per-pulse scaling applied.  Add AmpSF
% field if not already included.
if ~strcmp(cphd_meta.Data.SignalArrayFormat,'CF8')&&~isfield(cphd_meta.PVP,'AmpSF')
    cphd_meta.PVP.AmpSF = struct('Offset',0,'Size',1,'Format','F8');
    cphd_meta.Data.NumBytesPVP = cphd_meta.Data.NumBytesPVP + 8;
end
% Used internally within toolbox, but not official CRSD/CPHD field
if isfield(cphd_meta,'extra')
    cphd_meta = rmfield(cphd_meta,'extra');
end

%% Write file header and XML metadata
wb_handle = waitbar(0,['Converting to ' format]);
cphd_fid = fopen(output_cphd_file,'w','b'); % All CPHD must be big-endian
xml_string = sicdstruct2xml(cphd_meta, 'file_type', format, 'inc_class_attributes', false);
fh_struct = write_cphd_fileheader(cphd_fid, cphd_meta, format, numel(xml_string));
fwrite(cphd_fid,sprintf('\f\n'),'char');
fwrite(cphd_fid,sprintf('%s\f\n',xml_string),'char');

% Vector-based metadata is actually positioned before CPHD data in the
% file, but we have to write CPHD data first since we need to compute AmpSF
% before we can write vector-based metadata-- unless we want to read the
% entire input data twice.  However, we still have to allocate space in the
% file for vector-based metadata, since the fseek to SIGNAL_BLOCK_BYTE_OFFSET line
% below won't work to a location past the current end of the file.
fwrite(cphd_fid, zeros(fh_struct.PVP_BLOCK_SIZE + 2,1), 'uint8'); % Allocate space for vector-based metadata

%% Write CPHD data
% Iterate through pulses
waitbar(0,wb_handle,['Converting to ' format ' X (wideband)']);
fseek(cphd_fid,fh_struct.SIGNAL_BLOCK_BYTE_OFFSET,'bof'); % We should already be here.  Just for clarity of code
AmpSF = ones(length(p1.Results.pulse_indices), length(p1.Results.channels));
for i = 1:length(p1.Results.channels)
    for j = 1:length(p1.Results.pulse_indices)
        current_pulse = read_function(p1.Results.pulse_indices(j),'all',p1.Results.channels(i));
        pulse_interleaved = zeros(2*length(current_pulse),1); % FWRITE doesn't
        pulse_interleaved(1:2:end)=real(current_pulse); % handle complex data.
        pulse_interleaved(2:2:end)=imag(current_pulse); % Interleave I and Q.
        % We ignore any AmpSF from the input file and calculate our own.
        if ~strcmp(matlab_type,'float32') % Scale by AmpSF
            AmpSF(j,i) = max(pulse_interleaved)/double(intmax(matlab_type));
            pulse_interleaved = pulse_interleaved/AmpSF(j,i);
        end
        fwrite(cphd_fid,pulse_interleaved,matlab_type);
        waitbar(j/double(cphd_meta.Data.Channel(i).NumVectors),wb_handle);
    end
end

%% Write vector-based metadata
waitbar(0,wb_handle,['Converting to ' format ' X (narrowband)']);
fseek(cphd_fid,fh_struct.PVP_BLOCK_BYTE_OFFSET,'bof');
for i = 1:length(p1.Results.channels)
    [wbdata, nbdata] = read_function(p1.Results.pulse_indices, [], p1.Results.channels(i));
    if isfield(cphd_meta.PVP,'AmpSF')
        nbdata.AmpSF = AmpSF(:,i);
    end
    write_cphd_vbmeta(cphd_fid, cphd_meta, nbdata);
end
fwrite(cphd_fid,sprintf('\f\n'),'char');

%% Clean up, close input and output files
ph_reader.close();
fclose(cphd_fid);
close(wb_handle);

end

function fileheader = write_cphd_fileheader(cphd_fid, cphd_meta, format, xml_block_size)
% Get data element size
switch cphd_meta.Data.SignalArrayFormat
    case 'CI2'
        datatype_bytes=1;
    case 'CI4'
        datatype_bytes=2;
    case 'CF8'
        datatype_bytes=4;
end
% We don't know the offsets until we know the size of the file header.
% However, we don't know the file header until we know the size of the
% offsets.  A quandary...
fileheader.XML_BLOCK_SIZE = xml_block_size;
fileheader.XML_BLOCK_BYTE_OFFSET = 0; % Don't know this until we know the length of the file header
fileheader.PVP_BLOCK_SIZE = 0;
fileheader.PVP_BLOCK_BYTE_OFFSET = 0; % Don't know this yet
fileheader.SIGNAL_BLOCK_SIZE = 0;
fileheader.SIGNAL_BLOCK_BYTE_OFFSET = 0; % Don't know this yet
for channel=1:numel(cphd_meta.Data.Channel)
    fileheader.PVP_BLOCK_SIZE = fileheader.PVP_BLOCK_SIZE + ...
        cphd_meta.Data.NumBytesPVP*cphd_meta.Data.Channel(channel).NumVectors;
    fileheader.SIGNAL_BLOCK_SIZE = fileheader.SIGNAL_BLOCK_SIZE + ...
        2 * datatype_bytes *cphd_meta.Data.Channel(channel).NumVectors * ...
        cphd_meta.Data.Channel(channel).NumSamples;
end
if isfield(cphd_meta,'CollectionID')
    if isfield(cphd_meta.CollectionID,'Classification')
        fileheader.CLASSIFICATION = cphd_meta.CollectionID.Classification;
    end
    if isfield(cphd_meta.CollectionID,'ReleaseInfo')
        fileheader.RELEASE_INFO = cphd_meta.CollectionID.ReleaseInfo;
    end
end

% Potentially measuring the length of the fileheader could change the
% offsets for the XML, VB metadata, and CPHD (which are all originally
% measured off of a file header length of 0.) In the pathological case, it
% could take 5 times through this loop for all the string length of the
% file header to percolate through to all of the offsets.
fileheader_string = create_fh_string(struct(), format);
while fileheader.XML_BLOCK_BYTE_OFFSET ~= (length(fileheader_string))
    fileheader.XML_BLOCK_BYTE_OFFSET = length(fileheader_string);
    fileheader.PVP_BLOCK_BYTE_OFFSET = fileheader.XML_BLOCK_BYTE_OFFSET + fileheader.XML_BLOCK_SIZE + 2;
    fileheader.SIGNAL_BLOCK_BYTE_OFFSET = fileheader.PVP_BLOCK_BYTE_OFFSET + fileheader.PVP_BLOCK_SIZE;
    fileheader_string = create_fh_string(fileheader, format);
end

% After we have finalized all offsets, write string to file
fwrite(cphd_fid,fileheader_string,'char');

    % Create file header string
    function fh_str = create_fh_string(file_header_struct, format)
        if strcmp(format,'CRSD')
            fh_str = sprintf('CRSD/1.0.0\n'); % File Type Header
        elseif strcmp(format,'CPHD')
            fh_str = sprintf('CPHD/1.1.0\n'); % File Type Header
        end
        fileheader_fieldnames = fieldnames(file_header_struct);
        for j = 1:length(fileheader_fieldnames)
            % Key-Value Pair List
            if isnumeric(file_header_struct.(fileheader_fieldnames{j}))
                value = num2str(file_header_struct.(fileheader_fieldnames{j}));
            else 
                value = file_header_struct.(fileheader_fieldnames{j});
            end
            fh_str = sprintf('%s%s := %s\n',fh_str, ...
                fileheader_fieldnames{j},value);
        end
        fh_str = sprintf('%s\f\n',fh_str); % Section terminator
        
    end
end

function write_cphd_vbmeta(cphd_fid, cphd_meta, nbdata)
    % Convert vector-based data from structure to array, so it can be
    % written in a single FWRITE call.
    vectorParametersFlatten = cphd_meta.PVP;  % Copy structure
    % Flatten CRSD structures
    if isfield(vectorParametersFlatten,'RcvAntenna')
        vectorParametersFlatten = rmfield(vectorParametersFlatten, 'RcvAntenna');
        vectorParametersFlatten = setstructfields(vectorParametersFlatten, cphd_meta.PVP.RcvAntenna);
    end
    if isfield(vectorParametersFlatten,'TxPulse')
        vectorParametersFlatten = rmfield(vectorParametersFlatten, 'TxPulse');
        vectorParametersFlatten = setstructfields(vectorParametersFlatten, cphd_meta.PVP.TxPulse);
    end
    if isfield(vectorParametersFlatten,'TxAntenna')
        vectorParametersFlatten = rmfield(vectorParametersFlatten, 'TxAntenna');
        vectorParametersFlatten = setstructfields(vectorParametersFlatten, cphd_meta.PVP.TxPulse.TxAntenna);
    end
    vectorParametersCell = fieldnames(vectorParametersFlatten);
    vb_array = zeros(cphd_meta.Data.NumBytesPVP/8,length(nbdata.TxTime));
    for i = 1:length(vectorParametersCell)
        vb_array(vectorParametersFlatten.(vectorParametersCell{i}).Offset+...
            (1:vectorParametersFlatten.(vectorParametersCell{i}).Size),:) = ...
            nbdata.(vectorParametersCell{i}).';
    end
    fwrite(cphd_fid, vb_array, 'double');
end

        
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////