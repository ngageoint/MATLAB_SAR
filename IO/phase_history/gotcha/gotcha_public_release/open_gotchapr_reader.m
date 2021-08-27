function readerobj = open_gotchapr_reader( filename )
%OPEN_GOTCHAPR_READER Reads in data from the AFRL 2D/3D SAR Volumetric
% public release dataset as if it were CPHD "X" file format.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Read in all data in given directory

% Any file in the directory will work.  This reader will use all files in
% this directory (the entire pass of data) regardless of which file was
% actually selected.
[path, name] = fileparts(filename);

% The LOAD functions below brings in a variable of the name 'data'.
% Because of the helper subfunction, this is a "static workspace" and
% variables cannot be added at runtime. Therefore we pre-allocate this
% variable name.
data = [];

% Its much actually faster to read all the data twice, so that we can first
% determine the total number of pulses in the pass and then preallocate
% arrays for all data.
total_num_pulses = 0;
for az = 1:360
    name(20:22) = sprintf('%03d',az);
    load(fullfile(path,name)); % Loads data into a structure named 'data'
    total_num_pulses = total_num_pulses + size(data.fp,2);
end

% Preallocate arrays for metadata and pulse data
% GOTCHA data is defined in scene-center oriented coordinate system.
% Technically CPHD should be in ECEF coords.  In the public release, the
% actually scene center coordinates were not given, so we put scene center
% at the north pole to make it "ECEF-like" so that at least IFPs computing
% the WGS84 tangent plane will calculate something reasonable.
b=6356752.314245179; % Semi-minor (polar) axis of WGS_84 model
SRP_ECEF = [0 0 b]; % North pole
all_nbdata.TxTime = zeros(total_num_pulses,1); % Not available
all_nbdata.TxPos = zeros(total_num_pulses,3);
all_nbdata.RcvTime = zeros(total_num_pulses,1); % Not available
all_nbdata.RcvPos = zeros(total_num_pulses,3);
all_nbdata.SRPPos = repmat(SRP_ECEF,total_num_pulses,1);
all_nbdata.SC0 = zeros(total_num_pulses,1);
all_nbdata.SCSS = zeros(total_num_pulses,1);
all_nbdata.FX1 = min(data.freq)*ones(total_num_pulses,1);
all_nbdata.FX2 = max(data.freq)*ones(total_num_pulses,1);
pulse_data = zeros(size(data.fp,1),total_num_pulses);

% This reader will read and hold all pulse data in memory upon opening.
% This is somewhat contrary to the general concept of the phase history
% reader object, which can generally selectively read pulse data from a
% file.  However, this dataset is small enough that we don't worry about
% that here. (An entire circle pass easily fits in memory and only takes
% about 1 second to read on the machine this function was developed on.)
last_pulse = 0;
for az = 1:360
    name(20:22) = sprintf('%03d',az);
    load(fullfile(path,name)); % Load data in variable named 'data'
    num_pulses = size(data.fp,2);
    all_nbdata.TxPos(last_pulse+(1:num_pulses),:) = [data.x' data.y' (data.z + b)'];
    all_nbdata.SC0(last_pulse+(1:num_pulses)) = data.freq(1);
    all_nbdata.SCSS(last_pulse+(1:num_pulses)) = mean(diff(data.freq));
    pulse_data(:,last_pulse+(1:num_pulses)) = data.fp;
    
    last_pulse = last_pulse + num_pulses;
end
all_nbdata.RcvPos = all_nbdata.TxPos; % Unique positions not given, so monostatic case assumed

%% Setup CPHD XML metadata.
% CollectionID
cphd_meta.CollectionID.CollectorName = 'GOTCHA';
cphd_meta.CollectionID.CoreName = 'GOTCHA_public_release_dataset';
cphd_meta.CollectionID.CollectType = 'MONOSTATIC';
cphd_meta.CollectionID.RadarMode.ModeType = 'SPOTLIGHT';
cphd_meta.CollectionID.Classification = 'UNCLASSIFIED';

% Data
cphd_meta.Data.SignalArrayFormat = 'CF8';
cphd_meta.Data.NumCPHDChannels = 1;
% cphd_meta.Data.NumBytesPVP = 120;  % Let this be derived
cphd_meta.Data.Channel.NumVectors = total_num_pulses;
cphd_meta.Data.Channel.NumSamples = size(pulse_data,1);
cphd_meta.Data.Channel.SignalArrayByteOffset = 0;
cphd_meta.Data.Channel.PVPArrayByteOffset = 0;

% Global
cphd_meta.Global.DomainType = 'FX';
cphd_meta.Global.SGN = -1;
cphd_meta.Global.Timeline.CollectStart = datenum(2006,1,1); % Actual date not given in public release data, but we know the year
% The following times are totally bogus.  No time info given in this dataset.
cphd_meta.Global.Timeline.TxTime1 = 0;
cphd_meta.Global.Timeline.TxTime2 = 1;
cphd_meta.Global.FxBand.FxMin=min(all_nbdata.SC0);
cphd_meta.Global.FxBand.FxMax=max(all_nbdata.SC0+(all_nbdata.SCSS*size(pulse_data,1)));

% Channel
cphd_meta.Channel.FXFixedCPHD = true;
cphd_meta.Channel.TOAFixedCPHD = true;
cphd_meta.Channel.SRPFixedCPHD = true;
cphd_meta.Channel.Parameters.FXFixed = true;
cphd_meta.Channel.Parameters.TOAFixed = true;
cphd_meta.Channel.Parameters.SRPFixed = true;
cphd_meta.Channel.Parameters.SignalNormal = true;
cphd_meta.Channel.Parameters.FxC = mean(all_nbdata.SC0) + ...
    (mean(all_nbdata.SCSS)*size(pulse_data,1)/2);
cphd_meta.Channel.Parameters.FxBW = mean(all_nbdata.SCSS)*size(pulse_data,1);
% Maximum swath width support by sampling rate
cphd_meta.Channel.Parameters.TOASaved = 1/mean(all_nbdata.SCSS);
    
% SRP
cphd_meta.ReferenceGeometry.SRP.ECF.X = SRP_ECEF(1);
cphd_meta.ReferenceGeometry.SRP.ECF.Y = SRP_ECEF(2);
cphd_meta.ReferenceGeometry.SRP.ECF.Z = SRP_ECEF(3);

% Antenna
% No antenna information given in this dataset, so this required structure
% cannot be filled out.


%% Setup reader object
readerobj.read_cphd=@read_data; 
readerobj.get_meta=@() cphd_meta;
readerobj.close=@() 1;

    %% Functino for READ_CPHD method
    function [wbvectors, nbdata] = read_data(pulse_indices, sample_indices, channels)
        % Parse input parameters
        if (nargin<1)||strcmpi(pulse_indices,'all')
            pulse_indices=1:cphd_meta.Data.Channel.NumVectors;
        end
        if (nargin<2)||strcmpi(sample_indices,'all')
            sample_indices=1:cphd_meta.Data.Channel.NumSamples;
        end
        % Only one channel per directory, so we ignore the channel
        % parameter.  We could in theory setup the reader so that each
        % polarization is a channel from a single collect, but we have not
        % done this here.
        
        % Not much to do here, since we have already read all the data from
        % file.  Just return requested data.
        wbvectors = pulse_data(sample_indices,pulse_indices);
        nb_fieldnames = fieldnames(all_nbdata);
        for i=1:length(nb_fieldnames)
            nbdata.(nb_fieldnames{i}) = all_nbdata.(nb_fieldnames{i})(pulse_indices,:);
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////