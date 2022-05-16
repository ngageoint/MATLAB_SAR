function [ sicd_meta ] = meta2sicd_palsar2led( native_meta )
%META2SICD_PALSAR2LED Converts ALOS PALSAR 2 leader file into a SICD-style metadata structure
%
% Takes as input a metadata structure from read_ceos_led_meta.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% CollectionInfo
if strcmp(native_meta.data.scene_id(1:5), 'ALOS2')
    sicd_meta.CollectionInfo.CollectorName = native_meta.data.scene_id(1:5);
elseif strcmp(native_meta.data.scene_id(1:5), 'STRIX')
    sicd_meta.CollectionInfo.CollectorName = native_meta.data.scene_id(1:6);
end
% Do we want to convert CoreName to NGA-style pattern?
sicd_meta.CollectionInfo.CoreName = strtrim(native_meta.data.scene_id);

%% GeoData
% TODO: EarthModel is actually a lie.  The data is actually provided in
% GRS80, not WGS84, but the two are pretty close, so we don't worry about
% the difference for now.
% sicd_meta.GeoData.EarthModel = native_meta.data.ellips % GRS80 not valid for SICD
sicd_meta.GeoData.EarthModel = 'WGS_84';
% These fields aren't filled out in level 1.1 data so this won't work:
% sicd_meta.GeoData.SCP.LLH.Lat = native_meta.data.geo_lat;
% sicd_meta.GeoData.SCP.LLH.Lon = native_meta.data.geo_long;
% sicd_meta.GeoData.SCP.LLH.HAE = native_meta.data.avg_terr;
% We will let GeoData.SCP get populated from the IMG metadata instead.

%% Grid
sicd_meta.Grid.ImagePlane = 'SLANT';
sicd_meta.Grid.Type = 'RGZERO';
sicd_meta.Grid.Row.SS = native_meta.data.pixel_spacing;
sicd_meta.Grid.Col.SS = native_meta.data.line_spacing;
sicd_meta.Grid.Row.Sgn = -1;
sicd_meta.Grid.Col.Sgn = -1;
fc = SPEED_OF_LIGHT() / native_meta.data.wavelength;
sicd_meta.Grid.Row.KCtr = 2 / native_meta.data.wavelength;
sicd_meta.Grid.Col.KCtr = 0;
sicd_meta.Grid.Row.ImpRespBW = 2 * native_meta.data.bw_rng * 1000 / SPEED_OF_LIGHT();
dop_bw = native_meta.data.bw_az; % Doppler bandwidth
ss_zd_s = 1000/native_meta.data.prf; % Zero doppler spacing in seconds.  Could also use line timing in img file.
sicd_meta.Grid.Col.ImpRespBW = dop_bw * ss_zd_s / sicd_meta.Grid.Col.SS; % Convert to azimuth spatial bandwidth (cycles per meter)
if isfinite(native_meta.data_qual.sr_res)
    sicd_meta.Grid.Row.ImpRespWid = native_meta.data_qual.sr_res;
end
if isfinite(native_meta.data_qual.az_res)
    sicd_meta.Grid.Col.ImpRespWid = native_meta.data_qual.az_res;
end
sicd_meta.Grid.Row.DeltaKCOAPoly = 0;
if strcmpi(strtrim(native_meta.data.wgt_az),'1')
    sicd_meta.Grid.Row.WgtType.WindowName = 'UNIFORM';
end
if strcmpi(strtrim(native_meta.data.wgt_rng),'1')
    sicd_meta.Grid.Col.WgtType.WindowName = 'UNIFORM';
end

%% Timeline
% We could do this here, but its probably better to do in
% meta2sicd_palsar2img.m.  These values should be very close.
% prf = native_meta.data.prf/1000;
% sicd_meta.Timeline.IPP.Set.IPPPoly=[0; prf];

%% RadarCollection
sicd_meta.RadarCollection.TxFrequency.Min = fc - native_meta.data.bw_rng * 1000/2;
sicd_meta.RadarCollection.TxFrequency.Max = fc + native_meta.data.bw_rng * 1000/2;
if strcmpi(native_meta.data.range_pulse_code,'LINEAR FM CHIRP')
    sicd_meta.RadarCollection.WaveForm.WFParameters.TxPulseLength = native_meta.data.pulse_width * 1e-6;
    sicd_meta.RadarCollection.WaveForm.WFParameters.TxRFBandwidth = native_meta.data.bw_rng * 1000;
    sicd_meta.RadarCollection.WaveForm.WFParameters.TxFreqStart = sicd_meta.RadarCollection.TxFrequency.Min;
    sicd_meta.RadarCollection.WaveForm.WFParameters.TxFMRate = native_meta.data.range_pulse_amp_coef2;
    sicd_meta.RadarCollection.WaveForm.WFParameters.RcvDemodType = 'CHIRP';
    % sicd_meta.RadarCollection.WaveForm.WFParameters.RcvWindowLength = ???;
    sicd_meta.RadarCollection.WaveForm.WFParameters.ADCSampleRate = native_meta.data.sampling_rate;
end

%% ImageFormation
% TxFrequencyProc is not given explicitly in native_meta, so it will be
% populated by default in derived_sicd_fields()
% NumChanProc here does not consider the "channels" from single-beam vs
% dual-beam, as the term "channel" is used in native_meta.data.sar_chan.
sicd_meta.ImageFormation.RcvChanProc.NumChanProc = 1;
sicd_meta.ImageFormation.ImageFormAlgo = 'RMA'; % We assume
% There is really very little information available on the ALOS PALSAR
% image formation processor, so these are just guesses.
sicd_meta.ImageFormation.STBeamComp = 'NO';
sicd_meta.ImageFormation.ImageBeamComp = 'NO';
if strcmpi(strtrim(native_meta.data.autofocus_flg), 'YES')
    sicd_meta.ImageFormation.AzAutofocus = 'GLOBAL';  % No idea whether this is GLOBAL or SV
else
    sicd_meta.ImageFormation.AzAutofocus = 'NO';
end
sicd_meta.ImageFormation.RgAutofocus = 'NO';

%% RMA
sicd_meta.RMA.RMAlgoType='OMEGA_K';
sicd_meta.RMA.ImageType='INCA';
sicd_meta.RMA.INCA.FreqZero=fc;

%% Radiometric
% Some validation of this would be nice, but we believe this is the correct
% math based on the PALSAR documentation
sicd_meta.Radiometric.SigmaZeroSFPoly = 10^((native_meta.rad.cal_factor-32)/10);
% derived_sicd_fields should be able to derive RCSSFPoly, BetaZeroSFPoly,
% and GammaZeroSFPoly from this later.

%% ErrorStatistics
if all(isfinite([native_meta.pos.rad_pos_err, native_meta.pos.at_pos_err, ...
        native_meta.pos.ct_pos_err, native_meta.pos.rad_vel_err, ...
        native_meta.pos.at_vel_err, native_meta.pos.ct_vel_err]))
    sicd_meta.ErrorStatistics.Components.PosVelErr.Frame = 'RIC_ECF';
    sicd_meta.ErrorStatistics.Components.PosVelErr.P1 = native_meta.pos.rad_pos_err;
    sicd_meta.ErrorStatistics.Components.PosVelErr.P2 = native_meta.pos.at_pos_err;
    sicd_meta.ErrorStatistics.Components.PosVelErr.P3 = native_meta.pos.ct_pos_err;
    sicd_meta.ErrorStatistics.Components.PosVelErr.V1 = native_meta.pos.rad_vel_err;
    sicd_meta.ErrorStatistics.Components.PosVelErr.V2 = native_meta.pos.at_vel_err;
    sicd_meta.ErrorStatistics.Components.PosVelErr.V3 = native_meta.pos.ct_vel_err;
    sicd_meta.ErrorStatistics.Components.RadarSensor.RangeBias = 0.01; % Don't know this.  Just put a small number.
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////