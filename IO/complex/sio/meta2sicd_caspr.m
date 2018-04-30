function [ output_meta ] = meta2sicd_caspr( input_meta )
%META2SICD_CASPR Converts CASPR header data from read_caspr_meta into a SICD-style structure
%
% There are a number of variants of CASPR.  This function attempts to
% flexibly handle several of them, but there are probably more that this
% code has not been tested on, and the code may require some small tweaking
% for those untested variants.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if isfield(input_meta,'SpeedofLightParameters')&&...
        strcmp(input_meta.SpeedofLightParameters.distanceunitsDUNITSfeetmeters,'feet')
    factor=1/FEET_TO_METERS; % Set meters-to-feet conversion constant
else
    factor=1;
end

% Must be a better (sensor-independent) way to determine PCS vs. ECEF?
pcs_coords = isfield(input_meta,'DCSMissionModeParameters');
if pcs_coords
    xaxis_orientation = input_meta.DCSMissionModeParameters.PCSXaxisorientationwithrespecttonorthdegcw;
end

fn_index = structfun(@(x) isfield(x,'ImageIdentifier'), input_meta);
if any(fn_index)
    in_fnames = fieldnames(input_meta);
    output_meta.CollectionInfo.CoreName = input_meta.(in_fnames{fn_index}).ImageIdentifier(1:16);
    output_meta.CollectionInfo.CountryCode = input_meta.(in_fnames{fn_index}).Segment1CountryCode;
    Year = input_meta.(in_fnames{fn_index}).ImagingCDPCollStartTimeYear;
    DOY = input_meta.(in_fnames{fn_index}).ImagingCDPCollStartTimeDay;
    Hour = input_meta.(in_fnames{fn_index}).ImagingCDPCollStartTimeHours;
    Min = input_meta.(in_fnames{fn_index}).ImagingCDPCollStartTimeMinutes;
    Sec = input_meta.(in_fnames{fn_index}).ImagingCDPCollStartTimeSeconds;
    output_meta.Timeline.CollectStart = datenum([Year 0 DOY Hour Min Sec]);
elseif isfield(input_meta, 'FilenamesAndPathnames')
    output_meta.CollectionInfo.CoreName = input_meta.FilenamesAndPathnames.IFProotfilename;
elseif isfield(input_meta, 'FilenamesandPathnames')
    output_meta.CollectionInfo.CoreName = input_meta.FilenamesandPathnames.IFProotfilename;
end
switch input_meta.RADARParameters.collectionmodeframescanstripmap
    case 'frame'
        output_meta.CollectionInfo.RadarMode.ModeType = 'SPOTLIGHT';
    case 'scan'
        output_meta.CollectionInfo.RadarMode.ModeType = 'DYNAMIC STRIPMAP';
    case 'stripmap'
        output_meta.CollectionInfo.RadarMode.ModeType = 'STRIPMAP';
end
if isfield(input_meta, 'classification')
    output_meta.CollectionInfo.Classification = input_meta.classification;
end
output_meta.ImageCreation.DateTime=datenum(input_meta.CASPRCommonHeaderParameters.dateofprocessing);

output_meta.GeoData.EarthModel='WGS_84'; % Constant for all SICD
output_meta.GeoData.SCP.ECF.X=input_meta.ScenePositionParameters.xcomponentoutputreferencepointORPDUNITSPCS/factor;
output_meta.GeoData.SCP.ECF.Y=input_meta.ScenePositionParameters.ycomponentoutputreferencepointORPDUNITSPCS/factor;
output_meta.GeoData.SCP.ECF.Z=input_meta.ScenePositionParameters.zcomponentoutputreferencepointORPDUNITSPCS/factor;
if isfield(input_meta.InformationalParameters,'ECEFXpositionofprocessedscenecenterDUNITS')
    output_meta.GeoData.SCP.ECF.X=input_meta.InformationalParameters.ECEFXpositionofprocessedscenecenterDUNITS/factor;
    output_meta.GeoData.SCP.ECF.Y=input_meta.InformationalParameters.ECEFYpositionofprocessedscenecenterDUNITS/factor;
    output_meta.GeoData.SCP.ECF.Z=input_meta.InformationalParameters.ECEFZpositionofprocessedscenecenterDUNITS/factor;
elseif isfield(input_meta.InformationalParameters,'ECEFXpositionofprocessedscenecentermtr')
    output_meta.GeoData.SCP.ECF.X=input_meta.InformationalParameters.ECEFXpositionofprocessedscenecentermtr/factor;
    output_meta.GeoData.SCP.ECF.Y=input_meta.InformationalParameters.ECEFYpositionofprocessedscenecentermtr/factor;
    output_meta.GeoData.SCP.ECF.Z=input_meta.InformationalParameters.ECEFZpositionofprocessedscenecentermtr/factor;
end
scp_ecf = [output_meta.GeoData.SCP.ECF.X output_meta.GeoData.SCP.ECF.Y output_meta.GeoData.SCP.ECF.Z];
output_meta.Position.GRPPoly=output_meta.GeoData.SCP.ECF;
output_meta.GeoData.SCP.LLH.Lat=input_meta.ScenePositionParameters.latitudeofprocessedscenecenterdegNpositive;
output_meta.GeoData.SCP.LLH.Lon=input_meta.ScenePositionParameters.longitudeofprocessedscenecenterdegEpositive;
output_meta.GeoData.SCP.LLH.HAE=input_meta.ScenePositionParameters.altitudeofprocessedscenecenterDUNITS/factor;

output_meta.Grid.ImagePlane=upper(input_meta.ImageParameters.imageoutputplaneslantgroundfocusother);
output_meta.Grid.Type='RGAZIM';
output_meta.Grid.Row.UVectECF.X=-input_meta.SlantPlaneBasisVectors.xcomponentofXdirectionslantplanebasisvectorPCS;
output_meta.Grid.Row.UVectECF.Y=-input_meta.SlantPlaneBasisVectors.ycomponentofXdirectionslantplanebasisvectorPCS;
output_meta.Grid.Row.UVectECF.Z=-input_meta.SlantPlaneBasisVectors.zcomponentofXdirectionslantplanebasisvectorPCS;
if pcs_coords
    output_meta.Grid.Row.UVectECF = convert_from_pcs(output_meta.Grid.Row.UVectECF, scp_ecf, xaxis_orientation);
end
output_meta.Grid.Row.ImpRespWid=input_meta.ImageParameters.rangeresolutioninslantplanewweightingDUNITS/factor;
output_meta.Grid.Row.ImpRespBW=0.886/(input_meta.ImageParameters.rangeresolutioninslantplanewoweightingDUNITS/factor);
input_meta.FinalResamplerparameters.finalresamplerenabled = ...
   ~isfield(input_meta,'FinalResamplerparameters') || ...
   ~isfield(input_meta.FinalResamplerparameters,'finalresamplerenabled') || ...
   strcmpi(input_meta.FinalResamplerparameters.finalresamplerenabled,'true');
if input_meta.FinalResamplerparameters.finalresamplerenabled
   output_meta.Grid.Row.SS=output_meta.Grid.Row.ImpRespWid/...
      input_meta.ImageParameters.finalimagesampledensityrangepixelsIPR;
else
   output_meta.Grid.Row.SS=output_meta.Grid.Row.ImpRespWid/...
      input_meta.ImageParameters.compressedimagesampledensityrangepixelsIPR;
end
fillrng=1/(output_meta.Grid.Row.ImpRespBW*output_meta.Grid.Row.SS);
output_meta.Grid.Row.DeltaK1=-1/(output_meta.Grid.Row.SS*(fillrng)*2);
output_meta.Grid.Row.DeltaK2=-output_meta.Grid.Row.DeltaK1;
if isfield(input_meta.WeightingParameters,'weightingtypeuniformtaylorsvaiiq2dwalter')
   output_meta.Grid.Row.WgtType.WindowName=upper(input_meta.WeightingParameters.weightingtypeuniformtaylorsvaiiq2dwalter);
elseif isfield(input_meta.WeightingParameters,'weightingtypeuniformtaylorsvaiiq2dwaltersegsva')
   output_meta.Grid.Row.WgtType.WindowName=upper(input_meta.WeightingParameters.weightingtypeuniformtaylorsvaiiq2dwaltersegsva);
end
output_meta.Grid.Col.UVectECF.X=-input_meta.SlantPlaneBasisVectors.xcomponentofYdirectionslantplanebasisvectorPCS;
output_meta.Grid.Col.UVectECF.Y=-input_meta.SlantPlaneBasisVectors.ycomponentofYdirectionslantplanebasisvectorPCS;
output_meta.Grid.Col.UVectECF.Z=-input_meta.SlantPlaneBasisVectors.zcomponentofYdirectionslantplanebasisvectorPCS;
if pcs_coords
    output_meta.Grid.Col.UVectECF = convert_from_pcs(output_meta.Grid.Col.UVectECF, scp_ecf, xaxis_orientation);
end
output_meta.Grid.Col.ImpRespWid=input_meta.ImageParameters.azimuthresolutioninslantplanewweightingDUNITS/factor;
output_meta.Grid.Col.ImpRespBW=0.886/(input_meta.ImageParameters.azimuthresolutioninslantplanewoweightingDUNITS/factor);
if input_meta.FinalResamplerparameters.finalresamplerenabled
   output_meta.Grid.Col.SS=output_meta.Grid.Col.ImpRespWid/...
      input_meta.ImageParameters.finalimagesampledensityazimuthpixelsIPR;
else
   output_meta.Grid.Col.SS=output_meta.Grid.Col.ImpRespWid/...
      input_meta.ImageParameters.compressedimagesampledensityazimuthpixelsIPR;
end
fillaz=1/(output_meta.Grid.Col.ImpRespBW*output_meta.Grid.Col.SS);
output_meta.Grid.Col.DeltaK1=-1/(output_meta.Grid.Col.SS*(fillaz)*2);
output_meta.Grid.Col.DeltaK2=-output_meta.Grid.Col.DeltaK1;
if isfield(output_meta.Grid.Row,'WgtType')
    output_meta.Grid.Col.WgtType=output_meta.Grid.Row.WgtType;
end

times = [input_meta.ImageParameters.timeoffirstpulsecontributingtoimagesec...
    input_meta.ImageParameters.timeofcenterpulsecontributingtoimagesec...
    input_meta.ImageParameters.timeoflastpulsecontributingtoimagesec];
times = times - times(1);
xpos = [input_meta.InformationalParameters.sensorxpositionataperturestartDUNITSPCS/factor...
    input_meta.InformationalParameters.sensorxpositionataperturecenterDUNITSPCS/factor...
    input_meta.InformationalParameters.sensorxpositionatapertureendDUNITSPCS/factor];
if isfield(input_meta.InformationalParameters,'sensorypositionatapertursestartDUNITSPCS')
    % Typo found in some CASPR files
    input_meta.InformationalParameters.sensorypositionataperturestartDUNITSPCS = ...
        input_meta.InformationalParameters.sensorypositionatapertursestartDUNITSPCS;
end
ypos = [input_meta.InformationalParameters.sensorypositionataperturestartDUNITSPCS/factor...
    input_meta.InformationalParameters.sensorypositionataperturecenterDUNITSPCS/factor...
    input_meta.InformationalParameters.sensorypositionatapertureendDUNITSPCS/factor];
zpos = [input_meta.InformationalParameters.sensorzpositionataperturestartDUNITSPCS/factor...
    input_meta.InformationalParameters.sensorzpositionataperturecenterDUNITSPCS/factor...
    input_meta.InformationalParameters.sensorzpositionatapertureendDUNITSPCS/factor];
if pcs_coords
    % Convert from PCS to ECEF
    ecf_pos = pcs_to_ecf([xpos; ypos; zpos],scp_ecf,xaxis_orientation,true);
    xpos = ecf_pos(1,:);
    ypos = ecf_pos(2,:);
    zpos = ecf_pos(3,:);
end
output_meta.Position.ARPPoly.X = fliplr(polyfit(times,xpos,2)).';
output_meta.Position.ARPPoly.Y = fliplr(polyfit(times,ypos,2)).';
output_meta.Position.ARPPoly.Z = fliplr(polyfit(times,zpos,2)).';

% output_meta.RadarCollection.RefFreqIndex=uint32(0); % Absence of this field means all frequencies are true values
fc=input_meta.RADARParameters.transmitradarcenterfrequencyHz;
output_meta.RadarCollection.Waveform.WFParameters.TxFMRate=input_meta.RADARParameters.chirprateHzsec;
output_meta.RadarCollection.Waveform.WFParameters.RcvFMRate = input_meta.RADARParameters.chirprateHzsec;
output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate=input_meta.RADARParameters.complexsamplingfrequencyHz;
output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth=...
    output_meta.RadarCollection.Waveform.WFParameters.TxFMRate*...
    input_meta.RADARParameters.ADsamplescollectedoverxmitpulsewidth/...
    output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate;
output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart=fc-...
    (output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth/2);
output_meta.RadarCollection.TxFrequency.Min=output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart;
output_meta.RadarCollection.TxFrequency.Max=...
    output_meta.RadarCollection.TxFrequency.Min+...
    output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth;

if isfield(input_meta.PolarInterpolationPIParameters,...
        'timeoffirstpulseactuallycontributingtopolargrid')
    timeoffirstpulse = input_meta.PolarInterpolationPIParameters.timeoffirstpulseactuallycontributingtopolargrid;
    timeoflastpulse = input_meta.PolarInterpolationPIParameters.timeoflastpulseactuallycontributingtopolargrid;
elseif isfield(input_meta.InformationalParameters,...
        'timeoffirstpulseactuallycontributingtopolargrid')
    timeoffirstpulse = input_meta.InformationalParameters.timeoffirstpulseactuallycontributingtopolargrid;
    timeoflastpulse = input_meta.InformationalParameters.timeoflastpulseactuallycontributingtopolargrid;
end

if isfield(input_meta.PolarInterpolationPIParameters,'PIgridcenterfrequencyrangeHz')
    output_meta.Grid.Row.KCtr = 2*input_meta.PolarInterpolationPIParameters.PIgridcenterfrequencyrangeHz/SPEED_OF_LIGHT;
end

output_meta.ImageFormation.ImageFormAlgo='PFA';
output_meta.ImageFormation.TStartProc = timeoffirstpulse -...
    input_meta.ImageParameters.timeoffirstpulsecontributingtoimagesec;
if output_meta.ImageFormation.TStartProc < 0
    output_meta.ImageFormation.TStartProc = 0;
end
output_meta.ImageFormation.TEndProc = timeoflastpulse -...
    input_meta.ImageParameters.timeoffirstpulsecontributingtoimagesec;
if isfield(input_meta,'IPSData')
    output_meta.ImageFormation.TxRcvPolarizationProc=upper([input_meta.IPSData.polarizationvvvhhhhv(1),':',input_meta.IPSData.polarizationvvvhhhhv(2)]);
end

output_meta.SCPCOA.SCPTime=input_meta.ImageParameters.timeofcenterpulsecontributingtoimagesec-...
    input_meta.ImageParameters.timeoffirstpulsecontributingtoimagesec;
if strcmp(input_meta.RADARParameters.collectionmodeframescanstripmap,'frame')
    output_meta.Grid.TimeCOAPoly=output_meta.SCPCOA.SCPTime;
end
output_meta.SCPCOA.ARPPos.X=xpos(2);
output_meta.SCPCOA.ARPPos.Y=ypos(2);
output_meta.SCPCOA.ARPPos.Z=zpos(2);
if isfield(input_meta.InformationalParameters,'sensorxvelocityataperturecenterDUNITSsPCS')
    output_meta.SCPCOA.ARPVel.X=input_meta.InformationalParameters.sensorxvelocityataperturecenterDUNITSsPCS/factor;
    output_meta.SCPCOA.ARPVel.Y=input_meta.InformationalParameters.sensoryvelocityataperturecenterDUNITSsPCS/factor;
    output_meta.SCPCOA.ARPVel.Z=input_meta.InformationalParameters.sensorzvelocityataperturecenterDUNITSsPCS/factor;
elseif isfield(input_meta.InformationalParameters,'sensorxvelocityataperturecenterDUNITSPCS')
    output_meta.SCPCOA.ARPVel.X=input_meta.InformationalParameters.sensorxvelocityataperturecenterDUNITSPCS/factor;
    output_meta.SCPCOA.ARPVel.Y=input_meta.InformationalParameters.sensoryvelocityataperturecenterDUNITSPCS/factor;
    output_meta.SCPCOA.ARPVel.Z=input_meta.InformationalParameters.sensorzvelocityataperturecenterDUNITSPCS/factor;
end
if pcs_coords
    output_meta.SCPCOA.ARPVel = convert_from_pcs(output_meta.SCPCOA.ARPVel, scp_ecf, xaxis_orientation);
end
output_meta.SCPCOA.SideOfTrack=upper(input_meta.ScenePositionParameters.radarlookdirectionleftright(1));
output_meta.SCPCOA.SlantRange=input_meta.ImageParameters.slantrangeataperturecenterDUNITS/factor;
output_meta.SCPCOA.DopplerConeAng=input_meta.ImageParameters.coneangleataperturecenterdeg;
output_meta.SCPCOA.GrazeAng=input_meta.ImageParameters.grazingangleataperturecenterdeg;
output_meta.SCPCOA.IncidenceAng=90-output_meta.SCPCOA.GrazeAng;
output_meta.SCPCOA.TwistAng=-input_meta.ImageParameters.slantplanetiltangledeg;
output_meta.SCPCOA.SlopeAng=input_meta.ImageParameters.slopeangleataperturecenterdeg;

if isfield(input_meta,'FocusPlaneFTPComponents')
    output_meta.PFA.FPN.X=input_meta.FocusPlaneFTPComponents.xcomponentoffocusplanenormalunitvectorPCS;
    output_meta.PFA.FPN.Y=input_meta.FocusPlaneFTPComponents.ycomponentoffocusplanenormalunitvectorPCS;
    output_meta.PFA.FPN.Z=input_meta.FocusPlaneFTPComponents.zcomponentoffocusplanenormalunitvectorPCS;
elseif isfield(input_meta,'FocusPlaneBasisVectors')
    output_meta.PFA.FPN.X=input_meta.FocusPlaneBasisVectors.xcomponentofZdirectionfocusplanebasisvectorPCS;
    output_meta.PFA.FPN.Y=input_meta.FocusPlaneBasisVectors.ycomponentofZdirectionfocusplanebasisvectorPCS;
    output_meta.PFA.FPN.Z=input_meta.FocusPlaneBasisVectors.zcomponentofZdirectionfocusplanebasisvectorPCS;
end
if pcs_coords
    output_meta.PFA.FPN = convert_from_pcs(output_meta.PFA.FPN, scp_ecf, xaxis_orientation);
end
output_meta.PFA.IPN.X=input_meta.SlantPlaneBasisVectors.xcomponentofZdirectionslantplanebasisvectorPCS;
output_meta.PFA.IPN.Y=input_meta.SlantPlaneBasisVectors.ycomponentofZdirectionslantplanebasisvectorPCS;
output_meta.PFA.IPN.Z=input_meta.SlantPlaneBasisVectors.zcomponentofZdirectionslantplanebasisvectorPCS;
if pcs_coords
    output_meta.PFA.IPN = convert_from_pcs(output_meta.PFA.IPN, scp_ecf, xaxis_orientation);
end

try
    output_meta.Timeline.CollectStart = datenum([...
        input_meta.CASPRCommonHeaderParameters.dateofcollection...
        input_meta.CASPRCommonHeaderParameters.timeofcollection],...
        'ddmmmyyyyHH:MM:SS');
end
output_meta.Timeline.CollectDuration = input_meta.ImageParameters.approximateintegrationtimesec;
output_meta.Timeline.IPP.Set.TStart=0;
output_meta.Timeline.IPP.Set.TEnd=output_meta.Timeline.CollectDuration;
output_meta.Timeline.IPP.Set.IPPStart=uint32(0);
output_meta.Timeline.IPP.Set.IPPEnd=uint32(input_meta.RADARParameters.effectivePRFHz*...
    output_meta.Timeline.CollectDuration);

if isfield(input_meta,'RCSCalibrationparameters') &&...
        isfield(input_meta.RCSCalibrationparameters,'RCScalibrationenabled') &&...
        strcmpi(input_meta.RCSCalibrationparameters.RCScalibrationenabled,'true')
    output_meta.Radiometric.RCSSFPoly = 1;
end

end

% Converts from PCS to ECEF coordinates using the SICD X/Y/Z structure
function ecef_struct = convert_from_pcs(pcs_struct, scp_ecf, xaxis_orientation)
    ecf_pos = pcs_to_ecf([pcs_struct.X pcs_struct.Y pcs_struct.Z],...
        scp_ecf, xaxis_orientation, false);
    ecef_struct = struct('X',ecf_pos(1),'Y',ecf_pos(2)','Z',ecf_pos(3));
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////