function [freq,mpsd_out] = MPSD(filename,PlotFlag,StartPulse,StopPulse,Stride,ProcessingPulses,Channel,wb)
%MPSD Computes the MPSD (Mean Power Spectral Density).  
%
% INPUTS:
%   filename  - optional : phase history filename or reader object from
%                          open_ph_reader
%               missing  : open file dialog box opened to get file name
%               empty    : open file dialog box opened to get file name
%   PlotFlag  - optional : plot flag (1=Plot MPSD, 0 = No Plot)
%               missing  : plot flag set to 0
%               empty    : plot flag set to 0
%   StartPulse- optional : Starting Pulse
%               missing  : Start Pulse set to 1
%               empty    : Start Pulse set to 1
%   StopPulse - optional : Stop Pulse
%               missing  : Stop Pulse set to NumPulses
%               empty    : Stop Pulse set to NumPulses
%   Stride    - optional : Stride length through pulse data
%               missing  : Stride length set to 1
%               empty    : Stride length set to 1
%   ProcessingPulses - optional : Number of pulses to process for each step
%               missing  : Stride length set to 100
%               empty    : Stride length set to 100
%   Channel   - optional : Channel to Process
%               missing  : Set to 1
%               empty    : Set to 1
%   wb        - optional : Waitbar display flag
%               missing  : Set to 1
%               empty    : Set to 1
% OUTPUTS:
%   freq     : array of frequencies for each sample (Hz).  Assumes this is
%              the same for all pulses.
%   mpsd_out : computed Mean Power Spectral Density (dB)
%   MPSD plot (optional) : plot of MPSD 
%
% VERSION:
%   1.0
%     - Tim Cox 20080826
%     - initial version
%   1.1
%     - Tim Cox 20090610
%     - added .stats file format 
%   1.2
%     - Wade Schwartzkopf 20110928
%     - Reorganized code and setup to use open_ph_reader for all phase
%     history file formats rather than multiple case statements.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('Channel','var')
    Channel = 1; % Assume first channel if multi-channel dataset
end

%% Ask for filename if not passed in
if ~exist( 'filename', 'var' )||isempty(filename)
    [ fname, pathname ] = uigetfile( sar_file_extensions('phd'), ...
        'Open phase history file' );
    if fname == 0
        fprintf(1, 'MPSD :: file selection canceled by user\n');
        return;
    end
    filename = fullfile(pathname, fname);
end

%% Open file and compute MPSD
if ischar(filename)
    reader_obj = open_ph_reader(filename);
else
    reader_obj = filename;
end
% Extract relevant metadata
meta = reader_obj.get_meta();
NumSamples = meta.Data.Channel(Channel).NumSamples;
IID = meta.CollectionID.CoreName;
if isfield(meta.Data,'NumCRSDChannels')  % Check raw first if available
    if ~meta.Channel.Parameters(Channel).DemodFixed
        error('MPSD:VARIABLE_DEMOD','Fixed demod required.');
    end
    [~, nbdata] = reader_obj.read_raw(1,1,Channel);
    NumChannels = meta.Data.NumCRSDChannels;
    sampling_rate =  meta.Channel.Parameters(Channel).Fs;
    chirp = (nbdata.FICRate==0);
    if chirp
        fc = meta.Channel.Parameters(Channel).F0Ref + nbdata.DFIC0;
        freq = linspace(fc-sampling_rate/2,fc+sampling_rate/2,meta.Data.Channel(Channel).NumSamples);
    else
        freq = meta.Channel.Parameters(Channel).F0Ref + nbdata.DFIC0 + ...
            ((0:(double(meta.Data.Channel(Channel).NumSamples)-1))' * nbdata.FICRate / sampling_rate);
    end
elseif isfield(meta.Data,'NumCPHDChannels')
    if ~meta.Channel.FXFixedCPHD
        error('MPSD:VARIABLE_FX','Fixed FX parameters are required.');
    end
    [~, nbdata] = reader_obj.read_cphd(1,1,Channel);
    NumChannels = meta.Data.NumCPHDChannels;
    freq = linspace(nbdata.SC0,nbdata.SC0+nbdata.SCSS*double(NumSamples-1),...
        meta.Data.Channel(Channel).NumSamples);
end

% Handle defaults for optional input arguments
if ~exist('StartPulse', 'var')||isempty(StartPulse)||(StartPulse <= 0)
    StartPulse = 1;
end
if ~exist('StopPulse', 'var')||isempty(StopPulse)||(StopPulse <= 0)
    StopPulse = meta.Data.Channel(Channel).NumVectors;
end
if ~exist('Stride', 'var')||isempty(Stride)
    Stride = 1;
end
if ~exist('ProcessingPulses', 'var')||isempty(ProcessingPulses)
    ProcessingPulses = 100;
end
if ~exist('Channel', 'var')||isempty(Channel)
    Channel = 1;
end
if (Channel > NumChannels)
    msgbox('Invalid Channel Selection');
    return;
end
if ~exist('wb', 'var')||isempty(wb)
    wb = 1;
end

% Read data and compute MPSD
TotPow = zeros(meta.Data.Channel(Channel).NumSamples,1);
pulse_indices = double(StartPulse):double(Stride):double(StopPulse);
if wb
    wbh = waitbar(0, 'Processing pulses');
end
PULSE_BLOCK_SIZE = ProcessingPulses;
num_pulse_sets = ceil(numel(pulse_indices)/PULSE_BLOCK_SIZE);
for i=1:num_pulse_sets
    if wb
        waitbar((i-1)/num_pulse_sets, wbh, sprintf('Processing Pulse Step %d of %d',i,num_pulse_sets));
    end
    current_pulse_block = pulse_indices((i-1)*PULSE_BLOCK_SIZE+1:min(i*PULSE_BLOCK_SIZE,end));
    if isfield(meta.Data,'NumCRSDChannels')  % Check raw first if available
        pulses = reader_obj.read_raw(current_pulse_block,1:NumSamples,Channel);
        pulses = bsxfun(@minus, pulses, mean(pulses)); % Remove bias per pulse
        if chirp
            pulses = fft(pulses);
            pulses(1,:) = pulses(2,:);
            pulses = fftshift(pulses);
            pulses(1,:) = pulses(2,:);
        else % stretch
            pulses = deskew_rvp(pulses, sampling_rate, nbdata.FICRate);
        end
    elseif isfield(meta.Data,'NumCPHDChannels')
        pulses = reader_obj.read_cphd(current_pulse_block,1:NumSamples,Channel);
        if ~strcmp(meta.Global.DomainType, 'FX')
            pulses = fft(pulses);
        end
    end
    TotPow = TotPow + sum(abs(pulses).^2,2);
end
if wb
    waitbar(1, wbh);
    close(wbh);
end
mpsd_out = TotPow./numel(pulse_indices);

if ischar(filename)
    reader_obj.close(); % Close reader
end

%% Convert MPSD to dB
mpsd_out(mpsd_out <= 0) = 1;
mpsd_out = 10.*log10(mpsd_out);
mpsd_out = mpsd_out - min(mpsd_out); % Scale such that min is at 0 dB

if isfield(meta.Data,'NumCRSDChannels') && chirp %get rid of any DC bias
    CenterSample = round(NumSamples/2);
    mpsd_out(CenterSample-5:CenterSample+5) = (mpsd_out(CenterSample-5)+mpsd_out(CenterSample+5))/2;
end

%% Plot MPSD if PlotFlag is set
if exist('PlotFlag','var')&&PlotFlag
    figure;
    plot(freq./10^6,mpsd_out);
    title(sprintf('MPSD for IID: %s',IID));
    xlabel('Frequency (MHz)');
    ylabel('Uncalibrated Power (dB)');
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////