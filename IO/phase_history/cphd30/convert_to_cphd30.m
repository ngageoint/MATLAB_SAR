function convert_to_cphd30( input_ph_filename, output_cphd_file, varargin )
%CONVERT_TO_CPHD30 Converts a phase history dataset to CPHD 3.0 file format
%
% CONVERT_TO_CPHD30(INPUT_PH_FILENAME, OUTPUT_CPHD_FILENAME, 'PropertyName', PropertyValue, ...)
%
%       Property name     Description
%       pulse_indices     Pulses from the input phase history dataset to be
%                         written to the output CPHD file.  Default is all
%                         pulses.
%       channels          Channel from the input phase history dataset to
%                         be written to the output CPHD file.  Must be a
%                         single channel.  Default is 1.
%       SampleType        Datatype to be written to the output CPHD file.
%                         Default is whatever the sample type was closest
%                         to in the input file.  Valid options in CPHDX
%                         version 3.0 are 'RE08I_IM08I', 'RE16I_IM16I', and
%                         'RE32F_IM32F'.
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
all_vectors = 1:cphd_meta.Data.Channel(1).NumVectors;
p1.addParamValue('pulse_indices', all_vectors, @(x) all(ismember(x,all_vectors)));
p1.addParamValue('channels', 1, @(x) isscalar(x)&&ismember(x,1:cphd_meta.Data.NumCPHDChannels));
p1.addParamValue('SampleType', cphd_meta.Data.SignalArrayFormat, @ischar);
% Switching domain type not currently supported, but we leave a placeholder for this option
% p1.addParamValue('DomainType', cphd_meta.Global.DomainType, ...
%     @(x) sum(strcmp(x, {'FX','TOA'})==1));
% Option to writes raw data instead of compensated into CPHD format
p1.FunctionName = mfilename;
p1.parse(varargin{:});
% If we are changing any data, update metadata
if ~any(strcmp('pulse_indices',p1.UsingDefaults))
    cphd_meta.Data.Channel(p1.Results.channels).NumVectors = length(p1.Results.pulse_indices);
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
% Integer types need appropriate per-pulse scaling applied.  Add AmpSF
% field if not already included.
include_AmpSF0 = ~strcmp(cphd_meta.Data.SignalArrayFormat,'CF8');

%% Write CPHD preamble
wb_handle=waitbar(0,'Converting to CPHD 3.0');
cphd_fid = fopen(output_cphd_file,'w'); % Open with system default endian
write_cphd_preamble(cphd_fid, meta2cphd30_cphdx(cphd_meta, p1.Results.channels), include_AmpSF0);

% Allocate space for narrowband data.  Narrowband data actually comes
% before wideband data in the file, but we must write it after wideband
% data in our code, since we might need the AmpSF0 values calculated during
% the wideband writing.
nb_byte_offset = ftell(cphd_fid);
bytes_per_vector = 112; % Channel (4), vector (4), SRPPos (24), TxTime (8), TxPos (24), RcvTime (8), RcvPos (24), Fx0 (8), Fx_SS(8)
if include_AmpSF0, bytes_per_vector = bytes_per_vector + 8; end; % AmpSF0
nb_size_bytes = bytes_per_vector*cphd_meta.Data.Channel(p1.Results.channels).NumVectors;
fwrite(cphd_fid, zeros(nb_size_bytes,1), 'uint8');

%% Write wideband pulse data
% Iterate through pulses
waitbar(0,wb_handle,'Converting to CPHD 3.0 (wideband)');
for i=1:length(p1.Results.pulse_indices)
    fwrite(cphd_fid,p1.Results.channels-1,'uint32');
    fwrite(cphd_fid,i-1,'uint32');
    current_pulse = ph_reader.read_cphd(p1.Results.pulse_indices(i),'all',p1.Results.channels);
    pulse_interleaved = zeros(2*length(current_pulse),1); % FWRITE doesn't
    pulse_interleaved(1:2:end)=real(current_pulse); % handle complex data.
    pulse_interleaved(2:2:end)=imag(current_pulse); % Interleave I and Q.
    if include_AmpSF0 % Scale by AmpSF0
        AmpSF(i) = max(pulse_interleaved)/double(intmax(matlab_type));
        pulse_interleaved = pulse_interleaved/AmpSF(i);
    end
    fwrite(cphd_fid,pulse_interleaved,matlab_type);
    waitbar(i/double(cphd_meta.Data.Channel(p1.Results.channels).NumVectors),wb_handle);
end

%% Write per pulse metadata
waitbar(0,wb_handle,'Converting to CPHD 3.0 (narrowband)');
fseek(cphd_fid, nb_byte_offset, 'bof'); % Go back to location for narrowband data
[wbdata, nbdata] = ph_reader.read_cphd(p1.Results.pulse_indices, [], p1.Results.channels);
if include_AmpSF0
    nbdata.AmpSF = AmpSF;
end
write_cphd_nb(cphd_fid, nbdata, p1.Results.channels, include_AmpSF0);

%% Close input and output files
ph_reader.close();
fclose(cphd_fid);
close(wb_handle);

end

%%WRITE_CPHD_PREAMBLE Write out CPHD preamble
% Assumes file is already open
function write_cphd_preamble( fid, input_meta, include_AmpSF0 )

fprintf(fid, 'BegPreamble\n');
% Put metadata here
cphd_fieldnames=fieldnames(input_meta);
for i=1:length(cphd_fieldnames)
    value = input_meta.(cphd_fieldnames{i});
    if isnumeric(value)
        if strcmp(cphd_fieldnames{i},'DateTime')
            value = datestr(value, 'yyyymmddHHMMSS');
        else
            value = num2str(value);
        end
    elseif islogical(value)
        if value
            value = 'Yes';
        else
            value = 'No';
        end
    end
    fprintf(fid, '%-29s %s\n', cphd_fieldnames{i}, value);
end
fprintf(fid, 'BegNBVector\n');
fprintf(fid, '\tChannelNumber\n');
fprintf(fid, '\tVectorNumber\n');
fprintf(fid, '\tSRP\n');
fprintf(fid, '\tTxPos\n');
fprintf(fid, '\tRcvPos\n');
fprintf(fid, '\tTxTime\n');
fprintf(fid, '\tRcvTime\n');
fprintf(fid, '\tFx0\n');
fprintf(fid, '\tFxStepSize\n');
if include_AmpSF0, fprintf(fid, '\tAmpSF0\n'); end;
fprintf(fid, 'EndNBVector\n');
fprintf(fid, 'BegWBVector\n');
fprintf(fid, '\tChannelNumber\n');
fprintf(fid, '\tVectorNumber\n');
fprintf(fid, '\tSampleData\n');
fprintf(fid, 'EndWBVector\n');
fprintf(fid, 'EndPreamble\n');
fwrite(fid, hex2dec('4D6F6A6F'), 'uint32');

end

%%WRITE_CPHD_NB Writes narrowband metadata for CPHD file.  Assumes
%non-interleaved format.
function write_cphd_nb( fid, nb_meta, channel, include_AmpSF0 )

nb_meta.VectorNumber = (1:length(nb_meta.TxTime)).'; % Renumber pulses
for i=1:length(nb_meta.TxTime)
    fwrite(fid,channel,'uint32');
    fwrite(fid,nb_meta.VectorNumber(i)-1,'uint32');
    fwrite(fid,nb_meta.SRPPos(i,:),'double');
    fwrite(fid,nb_meta.TxPos(i,:),'double');
    fwrite(fid,nb_meta.RcvPos(i,:),'double');
    fwrite(fid,nb_meta.TxTime(i),'double');
    fwrite(fid,nb_meta.RcvTime(i),'double');
    fwrite(fid,nb_meta.SC0(i),'double');
    fwrite(fid,nb_meta.SCSS(i),'double');
    if include_AmpSF0, fwrite(fid,nb_meta.AmpSF(i),'double'); end;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////