function [ output_meta ] = meta2sicd_cmetaa( cmetaa_struct )
%META2SICD_CMETAA Converts metadata stored in the CMETAA NITF TRE into a
% SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

output_meta.CollectionInfo.RadarMode.ModeID=cmetaa_struct.RD_MODE;

output_meta.GeoData.EarthModel='WGS_84'; % Constant for all SICD
if strcmp(cmetaa_struct.CG_MODEL,'ECEF')
    scp_ecf = [cmetaa_struct.CG_SCECN_X cmetaa_struct.CG_SCECN_Y cmetaa_struct.CG_SCECN_Z];
    output_meta.GeoData.SCP.ECF.X=scp_ecf(1);
    output_meta.GeoData.SCP.ECF.Y=scp_ecf(2);
    output_meta.GeoData.SCP.ECF.Z=scp_ecf(3);
    scp_llh=ecf_to_geodetic(scp_ecf);
    output_meta.GeoData.SCP.LLH.Lat=scp_llh(1);
    output_meta.GeoData.SCP.LLH.Lon=scp_llh(2);
    output_meta.GeoData.SCP.LLH.HAE=scp_llh(3);
elseif strcmp(cmetaa_struct.CG_MODEL,'WGS84')
    scp_llh = [cmetaa_struct.CG_SCECN_X cmetaa_struct.CG_SCECN_Y cmetaa_struct.CG_SCECN_Z];
    output_meta.GeoData.SCP.LLH.Lat=cmetaa_struct.CG_SCECN_X;
    output_meta.GeoData.SCP.LLH.Lon=cmetaa_struct.CG_SCECN_Y;
    output_meta.GeoData.SCP.LLH.HAE=cmetaa_struct.CG_SCECN_Z;
    scp_ecf=geodetic_to_ecf(scp_llh);
    output_meta.GeoData.SCP.ECF.X=scp_ecf(1);
    output_meta.GeoData.SCP.ECF.Y=scp_ecf(2);
    output_meta.GeoData.SCP.ECF.Z=scp_ecf(3);
end
% These fields describe spatial frequency, not SCP.
% output_meta.ImageData.SCPPixel.Row=cmetaa_struct.IF_DC_IS_COL;
% output_meta.ImageData.SCPPixel.Col=cmetaa_struct.IF_DC_IS_ROW;
% output_meta.Position.GRPPoly = output_meta.GeoData.SCP.ECF;
if strcmp(cmetaa_struct.CG_MAP_TYPE,'GEOD')
    % Probably not in the right order
    output_meta.GeoData.ImageCorners.ICP.FRFC.Lat = cmetaa_struct.CG_PATCH_LTCORUL;
    output_meta.GeoData.ImageCorners.ICP.FRFC.Lon = cmetaa_struct.CG_PATCH_LGCORUL;
    output_meta.GeoData.ImageCorners.ICP.FRLC.Lat = cmetaa_struct.CG_PATCH_LTCORUR;
    output_meta.GeoData.ImageCorners.ICP.FRLC.Lon = cmetaa_struct.CG_PATCH_LGCORUR;
    output_meta.GeoData.ImageCorners.ICP.LRLC.Lat = cmetaa_struct.CG_PATCH_LTCORLR;
    output_meta.GeoData.ImageCorners.ICP.LRLC.Lon = cmetaa_struct.CG_PATCH_LGCORLR;
    output_meta.GeoData.ImageCorners.ICP.LRFC.Lat = cmetaa_struct.CG_PATCH_LTCORLL;
    output_meta.GeoData.ImageCorners.ICP.LRFC.Lon = cmetaa_struct.CG_PATCH_LNGCOLL;
end

if upper(cmetaa_struct.CMPLX_SIGNAL_PLANE(1))=='S'
    output_meta.Grid.ImagePlane='SLANT';
elseif upper(cmetaa_struct.CMPLX_SIGNAL_PLANE(1))=='G'
    output_meta.Grid.ImagePlane='GROUND';
end
output_meta.Grid.Row.SS=cmetaa_struct.IF_RSS;
output_meta.Grid.Row.ImpRespWid=cmetaa_struct.IF_RGRES;
output_meta.Grid.Row.Sgn=-str2double([cmetaa_struct.IF_RFFTS '1']); % Convention is opposite of SICD
row_zp = cmetaa_struct.IF_RFFT_TOT/cmetaa_struct.IF_RFFT_SAMP; % Zeropad factor
output_meta.Grid.Row.ImpRespBW=1/(cmetaa_struct.IF_RSS*row_zp);
output_meta.Grid.Col.SS=cmetaa_struct.IF_AZSS;
output_meta.Grid.Col.ImpRespWid=cmetaa_struct.IF_AZRES;
output_meta.Grid.Col.Sgn=-str2double([cmetaa_struct.IF_AFFTS '1']); % Convention is opposite of SICD
col_zp = cmetaa_struct.IF_AZFFT_TOT/cmetaa_struct.IF_AZFFT_SAMP; % Zeropad factor
output_meta.Grid.Col.ImpRespBW=1/(cmetaa_struct.IF_AZSS*col_zp);
% If sample spacing were ground plane values:
% output_meta.Grid.Row.SS=cosd(cmetaa_struct.CG_GAAC)*cmetaa_struct.IF_RSS;
% output_meta.Grid.Row.ImpRespWid=cosd(cmetaa_struct.CG_GAAC)*cmetaa_struct.IF_RGRES;
% output_meta.Grid.Col.SS=cosd(cmetaa_struct.CG_TILT)*cmetaa_struct.IF_AZSS; % Tilt is optional so test for this
% output_meta.Grid.Col.ImpRespWid=cosd(cmetaa_struct.CG_TILT)*cmetaa_struct.IF_AZRES;
% DeltaK1/2 will be derived later so we don't need to do this here:
% output_meta.Grid.Row.DeltaK1=-output_meta.Grid.Row.ImpRespBW/2;
% output_meta.Grid.Row.DeltaK2=-output_meta.Grid.Row.DeltaK1;
% output_meta.Grid.Col.DeltaK1=-output_meta.Grid.Col.ImpRespBW/2;
% output_meta.Grid.Col.DeltaK2=-output_meta.Grid.Col.DeltaK1;
switch cmetaa_struct.CMPLX_WEIGHT
    case 'UWT'
        output_meta.Grid.Row.WgtType.WindowName = 'UNIFORM';
        output_meta.Grid.Col.WgtType.WindowName = 'UNIFORM';
    case 'HMW'
        output_meta.Grid.Row.WgtType.WindowName = 'HAMMING';
        output_meta.Grid.Col.WgtType.WindowName = 'HAMMING';
    case 'HNW'
        output_meta.Grid.Row.WgtType.WindowName = 'HANNING';
        output_meta.Grid.Col.WgtType.WindowName = 'HANNING';
    case 'TAY'
        output_meta.Grid.Row.WgtType.WindowName = 'TAYLOR';
        output_meta.Grid.Row.WgtType.Parameter{1}.name = 'SLL';
        output_meta.Grid.Row.WgtType.Parameter{1}.value = num2str(-cmetaa_struct.CMPLX_RNG_SLL);
        output_meta.Grid.Row.WgtType.Parameter{2}.name = 'NBAR';
        output_meta.Grid.Row.WgtType.Parameter{2}.value = num2str(cmetaa_struct.CMPLX_RNG_TAY_NBAR);
        output_meta.Grid.Col.WgtType.WindowName = 'TAYLOR';
        output_meta.Grid.Col.WgtType.Parameter{1}.name = 'SLL';
        output_meta.Grid.Col.WgtType.Parameter{1}.value = num2str(-cmetaa_struct.CMPLX_AZ_SLL);
        output_meta.Grid.Col.WgtType.Parameter{2}.name = 'NBAR';
        output_meta.Grid.Col.WgtType.Parameter{2}.value = num2str(cmetaa_struct.CMPLX_AZ_TAY_NBAR);
end
if isfield(output_meta.Grid.Row,'WgtType') && ...
        ~strcmp(output_meta.Grid.Row.WgtType.WindowName, 'UNIFORM')
    % sicdweight2fun() knows how to convert SICD WgtType into function
    try % Without signal processing toolbox, taylorwin won't work
        rowfun = sicdweight2fun(output_meta.Grid.Row);
        output_meta.Grid.Row.WgtFunct = rowfun(512);
        colfun = sicdweight2fun(output_meta.Grid.Col);
        output_meta.Grid.Col.WgtFunct = colfun(512);
    end
end

try % Might fail if date is misformed.
    output_meta.Timeline.CollectStart = ...
        datenum([cmetaa_struct.T_UTC_YYYYMMMDD cmetaa_struct.T_HHMMSSUTC],'yyyymmmddHHMMSS');
catch
    % Nothing to do.  Just continue without this field.
end
output_meta.Timeline.CollectDuration=cmetaa_struct.WF_CDP;
output_meta.Timeline.IPP.Set.TStart = 0;
output_meta.Timeline.IPP.Set.TEnd = cmetaa_struct.WF_CDP;
output_meta.Timeline.IPP.Set.IPPStart = 0;
output_meta.Timeline.IPP.Set.IPPEnd = floor(cmetaa_struct.WF_CDP*cmetaa_struct.WF_PRF);
output_meta.Timeline.IPP.Set.IPPPoly = [0; cmetaa_struct.WF_PRF];

% output_meta.RadarCollection.RefFreqIndex=uint32(0); % Absence of this field means all frequencies are true values
output_meta.RadarCollection.TxFrequency.Min=cmetaa_struct.WF_SRTFR;
output_meta.RadarCollection.TxFrequency.Max=cmetaa_struct.WF_ENDFR;
output_meta.RadarCollection.TxPolarization=upper(cmetaa_struct.POL_TR);
output_meta.RadarCollection.Waveform.WFParameters.TxPulseLength = cmetaa_struct.WF_WIDTH;
output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth = cmetaa_struct.WF_BW;
output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart = cmetaa_struct.WF_SRTFR;
output_meta.RadarCollection.Waveform.WFParameters.TxFMRate = cmetaa_struct.WF_CHRPRT * 1e12;
output_meta.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization=...
    [upper(cmetaa_struct.POL_TR) ':' upper(cmetaa_struct.POL_RE)];

output_meta.SCPCOA.SCPTime = cmetaa_struct.WF_CDP/2; % No precise info given in complex NITF, so we use rough estimate
if strcmp(cmetaa_struct.CG_MODEL,'ECEF')
    output_meta.SCPCOA.ARPPos.X=cmetaa_struct.CG_APCEN_X;
    output_meta.SCPCOA.ARPPos.Y=cmetaa_struct.CG_APCEN_Y;
    output_meta.SCPCOA.ARPPos.Z=cmetaa_struct.CG_APCEN_Z;
elseif strcmp(cmetaa_struct.CG_MODEL,'WGS84')
    xyz=geodetic_to_ecf([cmetaa_struct.CG_APCEN_X cmetaa_struct.CG_APCEN_Y cmetaa_struct.CG_APCEN_Z]);
    output_meta.SCPCOA.ARPPos.X=xyz(1);
    output_meta.SCPCOA.ARPPos.Y=xyz(2);
    output_meta.SCPCOA.ARPPos.Z=xyz(3);
end
% Position.ARPPoly can be populated later by derived_sicd_fields()
if isfinite(cmetaa_struct.CG_SNVEL_X)&&cmetaa_struct.CG_SNVEL_X~=0 && ...
        isfinite(cmetaa_struct.CG_SNVEL_Y)&&cmetaa_struct.CG_SNVEL_Y~=0 && ...
        isfinite(cmetaa_struct.CG_SNVEL_Z)&&cmetaa_struct.CG_SNVEL_Z~=0
    output_meta.SCPCOA.ARPVel.X=cmetaa_struct.CG_SNVEL_X;
    output_meta.SCPCOA.ARPVel.Y=cmetaa_struct.CG_SNVEL_Y;
    output_meta.SCPCOA.ARPVel.Z=cmetaa_struct.CG_SNVEL_Z;
end
if isfinite(cmetaa_struct.CG_SNACC_X)&&cmetaa_struct.CG_SNACC_X~=0 && ...
        isfinite(cmetaa_struct.CG_SNACC_Y)&&cmetaa_struct.CG_SNACC_Y~=0 && ...
        isfinite(cmetaa_struct.CG_SNACC_Z)&&cmetaa_struct.CG_SNACC_Z~=0
    output_meta.SCPCOA.ARPAcc.X=cmetaa_struct.CG_SNACC_X;
    output_meta.SCPCOA.ARPAcc.Y=cmetaa_struct.CG_SNACC_Y;
    output_meta.SCPCOA.ARPAcc.Z=cmetaa_struct.CG_SNACC_Z;
end
if ~(isfield(output_meta,'SCPCOA')&&all(isfield(output_meta.SCPCOA,{'ARPPos','ARPVel'}))&&...
        isfield(output_meta,'GeoData')&&isfield(output_meta.GeoData,'SCP')&&...
        isfield(output_meta.GeoData.SCP,'ECF'))
    output_meta.SCPCOA.SideOfTrack=cmetaa_struct.CG_LD;
    output_meta.SCPCOA.SlantRange=cmetaa_struct.CG_SRAC;
    output_meta.SCPCOA.DopplerConeAng=cmetaa_struct.CG_CAAC;
    output_meta.SCPCOA.GrazeAng=cmetaa_struct.CG_GAAC;
    output_meta.SCPCOA.IncidenceAng=90-output_meta.SCPCOA.GrazeAng;
    if cmetaa_struct.CG_TILT>0 % Optional field
        output_meta.SCPCOA.TwistAng=cmetaa_struct.CG_TILT;
    end
    if cmetaa_struct.CG_SLOPE>0 % Optional field
        output_meta.SCPCOA.SlopeAng=cmetaa_struct.CG_SLOPE;
    end
else
    % Otherwise derived_sicd_fields can just compute these fields.
end

output_meta.ImageFormation.TxRcvPolarizationProc=...
    [upper(cmetaa_struct.POL_TR) ':' upper(cmetaa_struct.POL_RE)];
if strcmp(cmetaa_struct.IF_PROCESS,'PF')
    output_meta.ImageFormation.ImageFormAlgo='PFA';
    % Although not stated explicitly in the CMETAA spec, in practice it
    % appears that these unit vectors are generally populated in NED,
    % rather than ECF.
    fpn_ned = [cmetaa_struct.CG_FPNUV_X cmetaa_struct.CG_FPNUV_Y cmetaa_struct.CG_FPNUV_Z];
    fpn_ecf = ned_to_ecf(fpn_ned, scp_ecf, false);
    output_meta.PFA.FPN=cell2struct(num2cell(fpn_ecf),{'X','Y','Z'});
    ipn_ned = [cmetaa_struct.CG_IDPNUVX cmetaa_struct.CG_IDPNUVY cmetaa_struct.CG_IDPNUVZ];
    ipn_ecf = ned_to_ecf(ipn_ned, scp_ecf, false);
    output_meta.PFA.IPN=cell2struct(num2cell(ipn_ecf),{'X','Y','Z'});
elseif any(strcmp(cmetaa_struct.IF_PROCESS,{'RM','CD'}))
    output_meta.ImageFormation.ImageFormAlgo='RMA';
end
% CMETAA can't express all the ImageFormation description possible in SICD,
% so some of this is just a guess.  Normally when we aren't certain about a
% field, we would just leave it out, but these fields are all required
% fields in SICD, so the XML will not validate unless we put something
% here.  If these are known more precisely for a specific sensor, that
% adjustment should be put in the function that applies the sensor-specific
% portions of the NITF-SICD conversion.
% Since we aren't told any differently, we just show all time and frequency
% as having been processed.
output_meta.ImageFormation.TStartProc = 0;
output_meta.ImageFormation.TEndProc = cmetaa_struct.WF_CDP;
output_meta.ImageFormation.TxFrequencyProc.MinProc = cmetaa_struct.WF_SRTFR;
output_meta.ImageFormation.TxFrequencyProc.MaxProc = cmetaa_struct.WF_ENDFR;
output_meta.ImageFormation.STBeamComp = 'NO'; % Not sure if CMETAA can express this
if cmetaa_struct.IF_BEAM_COMP(1)=='Y'
    output_meta.ImageFormation.ImageBeamComp = 'SV';
else
    output_meta.ImageFormation.ImageBeamComp = 'NO';
end
if cmetaa_struct.AF_TYPE{1} =='N'
    output_meta.ImageFormation.AzAutofocus = 'NO';
else
    % No way to distinguish between SV and GLOBAL in CMETAA, so we guess.
    output_meta.ImageFormation.AzAutofocus = 'SV';
end
output_meta.ImageFormation.RgAutofocus = 'NO'; % Not sure if CMETAA can express this

% TODO: Should IF_RANGE_DATA swap Grid.Row/Col and
% ImageData.SCPPixel.Row/Col?  Probably.

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////