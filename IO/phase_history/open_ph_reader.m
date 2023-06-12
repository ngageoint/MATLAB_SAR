function [ readerobj ] = open_ph_reader( filename, varargin )
%OPEN_PH_READER Open a generic reader for multiple formats of phase history data
%
% Returns a reader object structure with the following methods:
%
% [wbvectors, per_pulse_metadata] = readerobj.read_cphd(pulse_indices, sample_indices, channels)
%      If file format holds raw phase history, the pulse data in the file
%      is read and then the correct processing is done so that motion
%      compensated pulses and metadata are returned.  If the file format
%      holds motion compensated phase history, this data is returned
%      directly.  Either way, this function always returns data at a
%      consistent stage in the processing chain.
%      PULSE_INDICES: Vector with indices to all desired pulses.  Can
%         also use the string 'all' to specify all pulses.  Default is all
%         pulses.
%      SAMPLES_INDICES: Vector with indices to all desired pulses.  Can
%         also use the string 'all' to specify all samples.  Default is
%         all samples.
%      CHANNELS: Vector with indices to all desired channels.  Can
%         also use the string 'all' to specify all pulses.  Default is 1.
%      WBVECTORS: The returned phase history (in 'wbvectors') is a matrix
%         of requested pulses wherein each pulse is a column with the first
%         requested pulse in column 1, second requested pulse in column 2,
%         etc.
%                      pulses ->   
%         samples |  | 1 2             |
%                 v  | 1 2 ...         |
%                    | 1 2             |
%      PER_PULSE_METADATA: MATLAB structure contain metadata for each
%         individual pulse.
%
%      Examples:
%         Read all data from all pulses, samples, and channels
%         [wbvectors, nbdata] = readerobj.read_cphd(filename,'all','all','all');
%
%         Read every other pulse and samples 1-1000 from channel number 4
%         [wbvectors, nbdata] = readerobj.read_cphd(filename,...
%                                                   first_pulse:2:last_pulse,
%                                                   1:1000,
%                                                   4);
%
%         Read only the narrowband data.  Will read preamble and all per-pulse
%         metadata, but WBVECTORS will be empty since no samples were requested.
%         [wbvectors, nbdata] = readerobj.read_cphd(filename,'all',[]); 
%
% [wbvectors, per_pulse_metadata] = readerobj.read_raw(pulse_indices, sample_indices, channels)
%      Reads raw pulse data without any phase stabilization/motion
%      compensation.  Not available for file formats that area already
%      motion compensated.  Input and output arguments are the same as for
%      the READ_CPHD method.
% readerobj.get_meta() % Returns a MATLAB structure with the CPHD
%         "preamble" metadata
% readerobj.close() % Should always be run after all reading is finished
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Determine format and open
format_string=guess_ph_format(filename);
if isempty(format_string)
    error('OPEN_PH_READER:UnrecognizedFormat','Unable to determine format of phase history file.');
end
readerobj = eval(['open_' format_string '_reader(''' filename ''', varargin{:});']);

% Define the get_nbdata for each reader object
if exist(ro,'read_cphd')
    readerobj.get_nbdata = ro.read_cphd('all', []);
elseif exist(ro,'read_raw')
    readerobj.get_nbdata = ro.read_raw('all', []);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////