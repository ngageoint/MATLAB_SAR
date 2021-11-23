function [ meta, symmetry ] = meta2sicd_tiffcapella( filename )
%META2SICD_TIFFCAPELLA Compute SICD metadata for Capella TIFF
%
% Written by: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

SECONDS_IN_A_DAY = 24*60*60;

% Examine TIFF
tags = read_tiff_tags(filename);
% MATLAB can't read compressed, blocked, complex TIFF, so user must
% translate first, if provided this type of data
if tags{1}.Compression~=1
    error('META2SICD_TIFFCAPELLA:COMPRESSED_TIFF', ['Unable to parse '...
        'compressed, blocked, and complex TIFF.  Please first '...
        'translate via ''gdal_translate -co TILED=no [capella tiff] '...
        '[uncompressed output]''.']);
end

% Initialize outputs
meta=struct();
meta.native.tiff = tags{1};
meta.native.tiff.ImageDescription = jsondecode(deblank(meta.native.tiff.ImageDescription));
symmetry=[0 strcmpi(meta.native.tiff.ImageDescription.collect.radar.pointing,'left') 0];

%% CollectionInfo
meta.CollectionInfo.CollectorName = meta.native.tiff.ImageDescription.collect.platform;
[startTime, startTimeFrac] = datenum_w_frac(...
    meta.native.tiff.ImageDescription.collect.start_timestamp);
% meta.CollectionInfo.CoreName = [... % NGA-style CoreName and prefix
%         upper(datestr(startTime,'ddmmmyy')) ...
%         meta.CollectionInfo.CollectorName ...
%         upper(datestr(startTime,'HHMMSS'))];
meta.CollectionInfo.CoreName = meta.native.tiff.ImageDescription.collect.collect_id;  % Capella-style CoreName
meta.CollectionInfo.CollectType = 'MONOSTATIC';
switch meta.native.tiff.ImageDescription.collect.mode
    case 'spotlight'
        meta.CollectionInfo.RadarMode.ModeType = 'SPOTLIGHT';
    case 'stripmap'
        meta.CollectionInfo.RadarMode.ModeType = 'STRIPMAP';
    case 'sliding_spotlight'
        meta.CollectionInfo.RadarMode.ModeType = 'DYNAMIC STRIPMAP';
    otherwise
        error('META2SICD_TIFFCAPELLA:UNRECOGNIZED_MODE', ...
            ['Mode ' meta.native.tiff.ImageDescription.collect.mode ' not currently handled.']);
end
meta.CollectionInfo.Classification = 'UNCLASSIFIED';

%% ImageCreation
meta.ImageCreation.Application=meta.native.tiff.Software;
meta.ImageCreation.DateTime=datenum(meta.native.tiff.ImageDescription.processing_time,'yyyy-mm-ddTHH:MM:SS');
meta.ImageCreation.Profile='Prototype';

%% ImageData
% Assumes symmetry(3)==0
meta.ImageData.NumCols=uint32(meta.native.tiff.ImageLength);
meta.ImageData.NumRows=uint32(meta.native.tiff.ImageWidth);
meta.ImageData.FullImage=meta.ImageData;
meta.ImageData.FirstRow=uint32(0); meta.ImageData.FirstCol=uint32(0);
if strcmpi(meta.native.tiff.ImageDescription.collect.image.data_type, 'CInt16')
    meta.ImageData.PixelType='RE16I_IM16I';
end
% Source: From communication with Capella
meta.ImageData.SCPPixel.Col = floor(meta.ImageData.NumCols/2);
meta.ImageData.SCPPixel.Row = floor(meta.ImageData.NumRows/2);
if strcmpi(meta.native.tiff.ImageDescription.collect.radar.pointing,'left')
    meta.ImageData.SCPPixel.Col = meta.ImageData.NumCols - meta.ImageData.SCPPixel.Col - 1;
end
% Note that meta.native.tiff.ImageDescription.collect.image.center_pixel.center_time
% is merely (start_timestamp + stop_timestamp)/2.  It does not denote the
% zero Doppler time of what the center pixel should be.

%% GeoData
meta.GeoData.EarthModel='WGS_84';
meta.GeoData.SCP.ECF.X = meta.native.tiff.ImageDescription.collect.image.center_pixel.target_position(1);
meta.GeoData.SCP.ECF.Y = meta.native.tiff.ImageDescription.collect.image.center_pixel.target_position(2);
meta.GeoData.SCP.ECF.Z = meta.native.tiff.ImageDescription.collect.image.center_pixel.target_position(3);
pos_llh = ecf_to_geodetic([meta.GeoData.SCP.ECF.X meta.GeoData.SCP.ECF.Y meta.GeoData.SCP.ECF.Z]);
meta.GeoData.SCP.LLH.Lat = pos_llh(1);
meta.GeoData.SCP.LLH.Lon = pos_llh(2);
meta.GeoData.SCP.LLH.HAE = pos_llh(3);

%% State vectors
num_state_vectors=numel(meta.native.tiff.ImageDescription.collect.state.state_vectors);
state_vector_T  = zeros(1,num_state_vectors);
state_vector_T_frac  = zeros(1,num_state_vectors);
state_vector_X  = zeros(1,num_state_vectors);
state_vector_Y  = zeros(1,num_state_vectors);
state_vector_Z  = zeros(1,num_state_vectors);
for i=1:num_state_vectors
    [state_vector_T(i), state_vector_T_frac(i)] = datenum_w_frac(...
        meta.native.tiff.ImageDescription.collect.state.state_vectors(i).time);
    
    state_vector_X(i) = meta.native.tiff.ImageDescription.collect.state.state_vectors(i).position(1);
    state_vector_Y(i) = meta.native.tiff.ImageDescription.collect.state.state_vectors(i).position(2);
    state_vector_Z(i) = meta.native.tiff.ImageDescription.collect.state.state_vectors(i).position(3);
end
state_vector_T = round((state_vector_T-startTime)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (state_vector_T_frac-startTimeFrac); % Handle fractional seconds
% sv2poly.m shows ways to determine best polynomial order, but 5th is
% almost always best for spaceborne systems
polyorder = min(1, numel(state_vector_T) - 1);
P_x = polyfit(state_vector_T, state_vector_X, polyorder);
P_y = polyfit(state_vector_T, state_vector_Y, polyorder);
P_z = polyfit(state_vector_T, state_vector_Z, polyorder);
meta.Position.ARPPoly.X = P_x(end:-1:1).';
meta.Position.ARPPoly.Y = P_y(end:-1:1).';
meta.Position.ARPPoly.Z = P_z(end:-1:1).';

%% Grid
if strcmp(meta.native.tiff.ImageDescription.product_type,'SLC')
    % Capella generally uses this grid, even for backprojected data
    meta.Grid.ImagePlane = 'SLANT';
    meta.Grid.Type = 'RGZERO';
end
% center_pixel.center time is actually zero Doppler time
[coaTime, coaTimeFrac] = datenum_w_frac(...
    meta.native.tiff.ImageDescription.collect.image.center_pixel.center_time);
meta.Grid.TimeCOAPoly = round((coaTime-startTime)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (coaTimeFrac-startTimeFrac); % Handle fractional seconds;
% Capella format and SICD label rows/columns differently
meta.Grid.Row.SS = meta.native.tiff.ImageDescription.collect.image.image_geometry.delta_range_sample;
% meta.native.tiff.ImageDescription.collect.image.pixel_spacing_column is ground spacing
meta.Grid.Col.SS = meta.native.tiff.ImageDescription.collect.image.pixel_spacing_row;
meta.Grid.Row.Sgn = -1;
meta.Grid.Col.Sgn = -1;
bw = meta.native.tiff.ImageDescription.collect.radar.time_varying_parameters.pulse_bandwidth;
if isfield(meta.native.tiff.ImageDescription.collect.image,'processed_range_bandwidth')
    meta.Grid.Row.ImpRespBW = 2*meta.native.tiff.ImageDescription.collect.image.processed_range_bandwidth/SPEED_OF_LIGHT;
else
    meta.Grid.Row.ImpRespBW = 2*bw/SPEED_OF_LIGHT;
end
meta.Grid.Row.ImpRespWid =  meta.native.tiff.ImageDescription.collect.image.range_resolution;
meta.Grid.Col.ImpRespWid =  meta.native.tiff.ImageDescription.collect.image.azimuth_resolution;
dop_bw = meta.native.tiff.ImageDescription.collect.image.processed_azimuth_bandwidth; % Doppler bandwidth
ss_zd_s = meta.native.tiff.ImageDescription.collect.image.image_geometry.delta_line_time;
meta.Grid.Col.ImpRespBW = dop_bw*abs(ss_zd_s)/meta.Grid.Col.SS; % Convert to azimuth spatial bandwidth (cycles per meter);
fc = meta.native.tiff.ImageDescription.collect.radar.center_frequency; % Center frequency
meta.Grid.Row.KCtr = 2*fc/SPEED_OF_LIGHT;
meta.Grid.Col.KCtr = 0;
meta.Grid.Row.DeltaK1 = -meta.Grid.Row.ImpRespBW/2;
meta.Grid.Row.DeltaK2 = -meta.Grid.Row.DeltaK1;
meta.Grid.Row.DeltaKCOAPoly = 0;
meta.Grid.Col.DeltaK1 = -meta.Grid.Col.ImpRespBW/2;
meta.Grid.Col.DeltaK2 = -meta.Grid.Col.DeltaK1;
meta.Grid.Row.WgtType.WindowName = ...
    meta.native.tiff.ImageDescription.collect.image.range_window.name;
par_names = fieldnames(meta.native.tiff.ImageDescription.collect.image.range_window.parameters);
if ~isempty(par_names)
    meta.Grid.Row.WgtType.Parameter.name = par_names{1};
    meta.Grid.Row.WgtType.Parameter.value =  num2str( ...
        meta.native.tiff.ImageDescription.collect.image.range_window.parameters.(...
        meta.Grid.Row.WgtType.Parameter.name));
end
if strcmpi(meta.Grid.Row.WgtType.WindowName,'avci-nacaroglu')
    meta.Grid.Row.WgtFunct = avci_nacaroglu_window(32,...
        meta.native.tiff.ImageDescription.collect.image.range_window.parameters.alpha);
end
if strcmpi(meta.native.tiff.ImageDescription.collect.image.azimuth_window.name,'antenna-taper')
    meta.Grid.Col.WgtType.WindowName = 'UNIFORM';
    % Antenna pattern has not been removed and no additional weighting was
    % applied.  This is further indicated below in the ImageFormation
    % section by setting STBeamComp to 'NO'.
else
    meta.Grid.Col.WgtType.WindowName = ...
        meta.native.tiff.ImageDescription.collect.image.azimuth_window.name;
    par_names = fieldnames(meta.native.tiff.ImageDescription.collect.image.azimuth_window.parameters);
    if ~isempty(par_names)
        meta.Grid.Col.WgtType.Parameter.name = par_names{1};
        meta.Grid.Col.WgtType.Parameter.value =  num2str( ...
            meta.native.tiff.ImageDescription.collect.image.azimuth_window.parameters.(...
            meta.Grid.Col.WgtType.Parameter.name));
    end
end

%% Radar Collection
meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth = bw;
meta.RadarCollection.Waveform.WFParameters.TxPulseLength = ...
    meta.native.tiff.ImageDescription.collect.radar.time_varying_parameters.pulse_duration;
meta.RadarCollection.Waveform.WFParameters.RcvDemodType='CHIRP';
meta.RadarCollection.Waveform.WFParameters.ADCSampleRate = ...
    meta.native.tiff.ImageDescription.collect.radar.sampling_frequency;
meta.RadarCollection.Waveform.WFParameters.RcvFMRate = 0; % True for RcvDemodType='CHIRP'
meta.RadarCollection.TxFrequency.Min = fc-(bw/2); % fc calculated in Grid section
meta.RadarCollection.TxFrequency.Max = fc+(bw/2);
% Assumes pulse parts are exactly adjacent in bandwidth
meta.RadarCollection.Waveform.WFParameters.TxFreqStart = ...
    meta.RadarCollection.TxFrequency.Min;
% Polarization
meta.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization = ...
    [meta.native.tiff.ImageDescription.collect.radar.transmit_polarization ...
    ':' meta.native.tiff.ImageDescription.collect.radar.receive_polarization];
meta.RadarCollection.TxPolarization = ...
    meta.native.tiff.ImageDescription.collect.radar.transmit_polarization;

%% Timeline
% TODO: Include multiple PRFs
meta.Timeline.CollectStart = startTime + (startTimeFrac/SECONDS_IN_A_DAY);
prf = meta.native.tiff.ImageDescription.collect.radar.prf.prf;
[endTime, endTimeFrac] = datenum_w_frac(...
    meta.native.tiff.ImageDescription.collect.stop_timestamp);
meta.Timeline.CollectDuration = round((endTime-startTime)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (endTimeFrac-startTimeFrac); % Handle fractional seconds;
meta.Timeline.IPP.Set.TStart = 0;
meta.Timeline.IPP.Set.TEnd = 0; % Apply real value later.  Just a placeholder.
meta.Timeline.IPP.Set.IPPStart = uint32(0);
meta.Timeline.IPP.Set.IPPEnd = uint32(meta.Timeline.CollectDuration*prf);
meta.Timeline.IPP.Set.IPPPoly = [0; prf];
meta.Timeline.IPP.Set.TEnd = double(meta.Timeline.IPP.Set.IPPEnd)/prf;

%% Image Formation
meta.ImageFormation.RcvChanProc = struct('NumChanProc', uint32(1), ...
    'PRFScaleFactor', 1);
if strcmp(meta.native.tiff.ImageDescription.collect.image.algorithm,'backprojection')
    % meta.ImageFormation.ImageFormAlgo = meta.native.tiff.ImageDescription.collect.image.algorithm;
    meta.ImageFormation.ImageFormAlgo = 'OTHER';  % Stay within SICD spec as much as possible
end
meta.ImageFormation.TStartProc = 0;
meta.ImageFormation.TEndProc = meta.Timeline.CollectDuration;
meta.ImageFormation.TxFrequencyProc.MinProc = ...
    meta.RadarCollection.TxFrequency.Min;
meta.ImageFormation.TxFrequencyProc.MaxProc = ...
    meta.RadarCollection.TxFrequency.Max;
meta.ImageFormation.TxRcvPolarizationProc = ...
    meta.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization;
meta.ImageFormation.STBeamComp = 'NO';
meta.ImageFormation.ImageBeamComp = 'NO';
meta.ImageFormation.AzAutofocus = 'NO';
meta.ImageFormation.RgAutofocus = 'NO';
meta.ImageFormation.Processing.Type = 'Backprojection';
meta.ImageFormation.Processing.Applied = true;

%% RMA
% Capella data is often processed with backprojection, but is formed to a
% RGZERO grid, so we treat it as RMA/INCA.
meta.RMA.RMAlgoType = 'RG_DOP';
meta.RMA.ImageType = 'INCA';
near_range = meta.native.tiff.ImageDescription.collect.image.image_geometry.range_to_first_sample;
meta.RMA.INCA.R_CA_SCP = near_range + ...
    (double(meta.ImageData.SCPPixel.Row)*meta.Grid.Row.SS);
meta.RMA.INCA.FreqZero = fc;
[centerTime, centerTimeFrac] = datenum_w_frac(...
    meta.native.tiff.ImageDescription.collect.image.center_pixel.center_time);
[firstTime, firstTimeFrac] = datenum_w_frac(...
    meta.native.tiff.ImageDescription.collect.image.image_geometry.first_line_time);
zd_t_scp = round((centerTime-firstTime)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (centerTimeFrac-firstTimeFrac);
meta.RMA.INCA.TimeCAPoly = [zd_t_scp; ss_zd_s/meta.Grid.Col.SS];
pos_coefs = [P_x(:) P_y(:) P_z(:)];
% Velocity is derivate of position.
vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
vel_x = polyval(vel_coefs(:,1), meta.Grid.TimeCOAPoly(1,1));
vel_y = polyval(vel_coefs(:,2), meta.Grid.TimeCOAPoly(1,1));
vel_z = polyval(vel_coefs(:,3), meta.Grid.TimeCOAPoly(1,1));
vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
meta.RMA.INCA.DRateSFPoly = 1/(sqrt(vm_ca_sq(1))*meta.RMA.INCA.TimeCAPoly(2));
% TODO: Doppler centroid

%% Radiometric
if strcmp(meta.native.tiff.ImageDescription.collect.image.radiometry,'beta_nought')
    meta.Radiometric.BetaZeroSFPoly = meta.native.tiff.ImageDescription.collect.image.scale_factor^2;
else
    warning('META2SICD_TIFFCAPELLA:UNRECOGNIZED_CALIBRATION', ...
            ['Radiometry mode ' meta.native.tiff.ImageDescription.collect.image.radiometry ' not currently handled.']);
end

meta = derived_sicd_fields(meta);

if isfield(meta,'Radiometric') && isfield(meta.Radiometric,'SigmaZeroSFPoly')
    meta.Radiometric.NoiseLevel.NoiseLevelType = 'ABSOLUTE';
    % meta.Radiometric.NoiseLevel.NoisePoly = ...
    %     meta.native.tiff.ImageDescription.collect.image.nesz_peak - 10*log10(meta.Radiometric.SigmaZeroSFPoly(1));
    meta.Radiometric.NoiseLevel.NoisePoly = polyshift(...
        meta.native.tiff.ImageDescription.collect.image.nesz_polynomial.coefficients(end:-1:1),meta.RMA.INCA.R_CA_SCP);
    meta.Radiometric.NoiseLevel.NoisePoly(1) = meta.Radiometric.NoiseLevel.NoisePoly(1) -...
        10*log10(meta.Radiometric.SigmaZeroSFPoly(1));
end

end

% MATLAB's datenum function won't handle precise times down under a
% millisecond, because 1) It won't accept a format with more than 3 .FFF in
% the string description of the date format and 2) the resulting serial
% date number is stored in days from 00-JAN-0000 and just doesn't have
% enough bits to handle fractional seconds to the level we want.  Here we
% handle the fractional seconds separately so we can read date with the
% precision we need.
function [datenum_s, datenum_frac] = datenum_w_frac(datestring)
    datenum_s = datenum(datestring,'yyyy-mm-ddTHH:MM:SS');
    datenum_frac = str2double(regexp(datestring,'\.\d*','match'));
    if isnan(datenum_frac), datenum_frac = 0; end
end

% Avci-Nacaroglu Exponential
function out = avci_nacaroglu_window(n, alpha)
    if ~exist('alpha','var')
        alpha = 1.25;
    end
    n2 = floor(n/2);
    t = ((0:(n-1)) - n2)/n;
    out = exp(pi * alpha * (sqrt(1 - (2 * t).^2) - 1));
end


% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////