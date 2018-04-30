function convert_to_cphdx( input_ph_filename, output_cphd_file, varargin )
%CONVERT_TO_CPHDX Converts a phase history dataset to CPHD X file format
% (version 0.3)
%
% CONVERT_TO_CPHDX(INPUT_PH_FILENAME, OUTPUT_CPHD_FILENAME, 'PropertyName', PropertyValue, ...)
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
%                         to in the input file.  Valid options in CPHDX
%                         version 3.0 are 'RE08I_IM08I', 'RE16I_IM16I', and
%                         'RE32F_IM32F'.
%       cphdi             Write raw, rather than compensated, phase history
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
% Author: Wade Schwartzkopf, NGA/IDT
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
all_vectors = 1:cphd_meta.Data.ArraySize(1).NumVectors;
p1.addParamValue('pulse_indices', all_vectors, @(x) all(ismember(x,all_vectors)));
all_channels = 1:cphd_meta.Data.NumCPHDChannels;
p1.addParamValue('channels', all_channels, @(x) all(ismember(x,all_channels)));
p1.addParamValue('SampleType', cphd_meta.Data.SampleType, @ischar);
% Switching domain type not currently supported, but we leave a placeholder for this option
% p1.addParamValue('DomainType', cphd_meta.Global.DomainType, ...
%     @(x) sum(strcmp(x, {'FX','TOA'})==1));
% Option to writes raw data instead of compensated into CPHD format
p1.addParamValue('cphdi', false, @(x) isscalar(x)&&islogical(x)&&isfield(ph_reader,'read_raw'));
p1.FunctionName = mfilename;
p1.parse(varargin{:});
% If we are changing any data, update metadata
if ~any(strcmp('pulse_indices',p1.UsingDefaults))
    for i = 1:cphd_meta.Data.NumCPHDChannels
        cphd_meta.Data.ArraySize(i).NumVectors = length(p1.Results.pulse_indices);
    end
end
if ~any(strcmp('SampleType',p1.UsingDefaults))
    cphd_meta.Data.SampleType = p1.Results.SampleType;
end    
% Convert CPHD type string to MATLAB type string
switch cphd_meta.Data.SampleType
    case 'RE08I_IM08I'
        matlab_type = 'int8';
    case 'RE16I_IM16I'
        matlab_type = 'int16';
    case 'RE32F_IM32F'
        matlab_type = 'float32';
    otherwise
        error('CONVERT_TO_CPHDX:UNRECOGNIZED_DATATYPE','Unrecognized data type.');
end
% Integer types need appropriate per-pulse scaling applied.  Add AmpSF
% field if not already included.
if ~strcmp(cphd_meta.Data.SampleType,'RE32F_IM32F')&&~isfield(cphd_meta.VectorParameters,'AmpSF')
    cphd_meta.VectorParameters.AmpSF = 8;
    cphd_meta.Data.NumBytesVBP = cphd_meta.Data.NumBytesVBP + cphd_meta.VectorParameters.AmpSF;
end
% Are we reading raw or compensated data?
if p1.Results.cphdi
    read_function = ph_reader.read_raw;
else
    read_function = ph_reader.read_cphd;
end

%% Write file header and XML metadata
wb_handle = waitbar(0,'Converting to CPHD X');
cphd_fid = fopen(output_cphd_file,'w','b'); % All CPHDX must be big-endian
xml_string = sicdstruct2xml(cphd_meta, 'file_type', 'CPHD', 'inc_class_attributes', false);
fh_struct = write_cphdx_fileheader(cphd_fid, cphd_meta, xml_string);
fwrite(cphd_fid,sprintf('\f\n'),'char');
fwrite(cphd_fid,sprintf('%s\f\n',xml_string),'char');

% Vector-based metadata is actually positioned before CPHD data in the
% file, but we have to write CPHD data first since we need to compute AmpSF
% before we can write vector-based metadata-- unless we want to read the
% entire input data twice.  However, we still have to allocate space in the
% file for vector-based metadata, since the fseek to CPHD_BYTE_OFFSET line
% below won't work to a location past the current end of the file.
fwrite(cphd_fid, zeros(fh_struct.VB_DATA_SIZE + 2,1), 'uint8'); % Allocate space for vector-based metadata

%% Write CPHD data
% Iterate through pulses
waitbar(0,wb_handle,'Converting to CPHD X (wideband)');
fseek(cphd_fid,fh_struct.CPHD_BYTE_OFFSET,'bof'); % We should already be here.  Just for clarity of code
AmpSF = zeros(length(p1.Results.pulse_indices), length(p1.Results.channels));
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
        elseif isfield(cphd_meta.VectorParameters,'AmpSF')
            AmpSF(j,i) = 1; % Float requires no scaling.
        end
        fwrite(cphd_fid,pulse_interleaved,matlab_type);
        waitbar(j/double(cphd_meta.Data.ArraySize(i).NumVectors),wb_handle);
    end
end

%% Write vector-based metadata
waitbar(0,wb_handle,'Converting to CPHD X (narrowband)');
fseek(cphd_fid,fh_struct.VB_BYTE_OFFSET,'bof');
for i = 1:length(p1.Results.channels)
    [wbdata, nbdata] = read_function(p1.Results.pulse_indices, [], p1.Results.channels(i));
    if isfield(cphd_meta.VectorParameters,'AmpSF')
        nbdata.AmpSF = AmpSF(:,i);
    end
    write_cphdx_vbmeta(cphd_fid, cphd_meta, nbdata);
end
fwrite(cphd_fid,sprintf('\f\n'),'char');

%% Clean up, close input and output files
ph_reader.close();
fclose(cphd_fid);
close(wb_handle);

end

function fileheader = write_cphdx_fileheader(cphd_fid, cphd_meta, xml_string)
% Get data element size
switch cphd_meta.Data.SampleType
    case 'RE08I_IM08I'
        datatype_bytes=1;
    case 'RE16I_IM16I'
        datatype_bytes=2;
    case 'RE32F_IM32F'
        datatype_bytes=4;
end
% We don't know the offsets until we know the size of the file header.
% However, we don't know the file header until we know the size of the
% offsets.  A quandary...
fileheader.XML_DATA_SIZE = length(xml_string);
fileheader.XML_BYTE_OFFSET = 0; % Don't know this until we know the length of the file header
fileheader.VB_DATA_SIZE = 0;
fileheader.VB_BYTE_OFFSET = 0; % Don't know this yet
fileheader.CPHD_DATA_SIZE = 0;
fileheader.CPHD_BYTE_OFFSET = 0; % Don't know this yet
for channel=1:cphd_meta.Data.NumCPHDChannels
    fileheader.VB_DATA_SIZE = fileheader.VB_DATA_SIZE + ...
        cphd_meta.Data.NumBytesVBP*cphd_meta.Data.ArraySize(channel).NumVectors;
    fileheader.CPHD_DATA_SIZE = fileheader.CPHD_DATA_SIZE + ...
        2 * datatype_bytes *cphd_meta.Data.ArraySize(channel).NumVectors * ...
        cphd_meta.Data.ArraySize(channel).NumSamples;
end
fileheader.CLASSIFICATION = cphd_meta.CollectionInfo.Classification;

% Potentially measuring the length of the fileheader could change the
% offsets for the XML, VB metadata, and CPHD (which are all originally
% measured off of a file header length of 0.) In the pathological case, it
% could take 5 times through this loop for all the string length of the
% file header to percolate through to all of the offsets.
fileheader_string = '';
while fileheader.XML_BYTE_OFFSET ~= (length(fileheader_string) + 2);
    fileheader.XML_BYTE_OFFSET = length(fileheader_string) + 2;
    fileheader.VB_BYTE_OFFSET = fileheader.XML_BYTE_OFFSET + fileheader.XML_DATA_SIZE + 2;
    fileheader.CPHD_BYTE_OFFSET = fileheader.VB_BYTE_OFFSET + fileheader.VB_DATA_SIZE + 2;
    fileheader_string = create_fh_string(fileheader);
end

% After we have finalized all offsets, write string to file
fwrite(cphd_fid,fileheader_string,'char');

    % Create file header string
    function fh_str = create_fh_string(file_header_struct)
        fh_str = sprintf('CPHD/0.3\n'); % File Type Header
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

function write_cphdx_vbmeta(cphd_fid, cphd_meta, nbdata)
    vectorParametersFlatten = cphd_meta.VectorParameters;
    if ~strcmp(cphd_meta.Data.SampleType,'RE32F_IM32F')
        vectorParametersFlatten.AmpSF = 8;
    end
    if strcmp(cphd_meta.Global.DomainType, 'FX')
        vectorParametersFlatten = rmfield(vectorParametersFlatten, 'FxParameters');
        vectorParametersFlatten = setstructfields(vectorParametersFlatten, cphd_meta.VectorParameters.FxParameters);
    elseif strcmp(cphd_meta.Global.DomainType, 'TOA')
        vectorParametersFlatten = rmfield(vectorParametersFlatten, 'TOAParameters');
        vectorParametersFlatten = setstructfields(vectorParametersFlatten, cphd_meta.VectorParameters.TOAParameters);
    else
        error('OPEN_CPHDX_READER:UNRECOGNIZED_DOMAIN_TYPE','Unrecognized domain type.');
    end
    vectorParametersCell = fieldnames(vectorParametersFlatten);

    % Convert vector-based data from structure to array, so it can be
    % written in a single FWRITE call.
    vb_array = zeros(cphd_meta.Data.NumBytesVBP/8,length(nbdata.TxTime));
    field_index = 1;
    for i = 1:length(vectorParametersCell)
        field_size = vectorParametersFlatten.(vectorParametersCell{i})/8;
        vb_array(field_index:(field_index+field_size-1),:) = nbdata.(vectorParametersCell{i}).';
        field_index = field_index + field_size;
    end
    fwrite(cphd_fid, vb_array, 'double');
end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////