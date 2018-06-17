function [ output_meta ] = meta2sicd_tsxxml( xml_domnode, geo_xml )
%META2SICD_TSX Converts TerraSAR-X XML description into a SICD-style metadata structure
%
% Takes as input a Document Object Model (DOM) node from the main TSX XML
% descriptor file and the geoReference XML file.
%
% SICD fields not currently computed:
%    ValidData (although computation of this is shown in cosar_valid_data.m)
% 
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup
SECONDS_IN_A_DAY = 24*60*60;
if ~exist('geo_xml','var')
    geo_xml = '';
end
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

%% CollectionInfo
output_meta.CollectionInfo.CollectorName=char(xp.evaluate(...
    'level1Product/productInfo/missionInfo/mission',xml_domnode));
output_meta.CollectionInfo.CoreName=char(xp.evaluate(...
    'level1Product/productInfo/sceneInfo/sceneID',xml_domnode));
output_meta.CollectionInfo.CollectType='MONOSTATIC';
output_meta.CollectionInfo.RadarMode.ModeID=char(xp.evaluate(...
    'level1Product/productInfo/acquisitionInfo/imagingMode',xml_domnode));
if strncmpi(output_meta.CollectionInfo.RadarMode.ModeID,'ST',2)
    output_meta.CollectionInfo.RadarMode.ModeType='SPOTLIGHT';
elseif any(strncmpi(output_meta.CollectionInfo.RadarMode.ModeID,{'SL','HS'},2)) % SL (spotlight), HS (high-resolution spotlight)
    output_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
else % SM (stripmap), SC (ScanSAR)
    output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
end
output_meta.CollectionInfo.Classification='UNCLASSIFIED';

%% ImageCreation
output_meta.ImageCreation.DateTime=datenum(char(xp.evaluate(...
    'level1Product/generalHeader/generationTime',...
    xml_domnode)),'yyyy-mm-ddTHH:MM:SS.FFF');
output_meta.ImageCreation.Profile='Prototype';

%% ImageData
% Rows/columns switched since data stored "sideways" from SICD standard
output_meta.ImageData.NumRows=uint32(str2double(xp.evaluate('level1Product/productInfo/imageDataInfo/imageRaster/numberOfColumns',xml_domnode)));
output_meta.ImageData.NumCols=uint32(str2double(xp.evaluate('level1Product/productInfo/imageDataInfo/imageRaster/numberOfRows',xml_domnode)));
output_meta.ImageData.FullImage=output_meta.ImageData;
output_meta.ImageData.FirstRow=uint32(0); output_meta.ImageData.FirstCol=uint32(0);
output_meta.ImageData.PixelType='RE16I_IM16I';
% A number of different options for selecting scene center point in TSX.
if ~isempty(geo_xml) % If geoReference file is available
    % We prefer to choose a central point from the geoReference geolocation
    % grid since those points have more precise timings.
    num_grid_az_pts = str2double(char(xp.evaluate('geoReference/geolocationGrid/numberOfGridPoints/azimuth',geo_xml)));
    num_grid_rg_pts = str2double(char(xp.evaluate('geoReference/geolocationGrid/numberOfGridPoints/range',geo_xml)));
    ctr_grid_az = floor(num_grid_az_pts/2)+1;
    ctr_grid_rg = ceil(num_grid_rg_pts/2);
    scp_geo_node_str = ['geoReference/geolocationGrid/gridPoint[@irg="' ...
        num2str(ctr_grid_rg) '" and @iaz="' num2str(ctr_grid_az) '"]'];
    output_meta.ImageData.SCPPixel.Row=uint32(str2double(char(xp.evaluate(...
        [scp_geo_node_str '/col'], geo_xml)))) - 1;
    output_meta.ImageData.SCPPixel.Col=uint32(str2double(char(xp.evaluate(...
        [scp_geo_node_str '/row'], geo_xml)))) - 1;
else  % Main TSX XML defines a center point, but its metadata not as precise
    output_meta.ImageData.SCPPixel.Row=uint32(str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCenterCoord/refColumn',xml_domnode))-1);
    output_meta.ImageData.SCPPixel.Col=uint32(str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCenterCoord/refRow',xml_domnode))-1);
end

%% GeoData
output_meta.GeoData.EarthModel='WGS_84';
% Initially, we just seed this with a rough value.  Later we will put in
% something more precise.
if ~isempty(geo_xml) % If geoReference file is available
    output_meta.GeoData.SCP.LLH.Lat=str2double(char(xp.evaluate(...
        [scp_geo_node_str '/lat'], geo_xml)));
    output_meta.GeoData.SCP.LLH.Lon=str2double(char(xp.evaluate(...
        [scp_geo_node_str '/lon'], geo_xml)));
    output_meta.GeoData.SCP.LLH.HAE=str2double(char(xp.evaluate(...
        [scp_geo_node_str '/height'], geo_xml)));
else
    output_meta.GeoData.SCP.LLH.Lat=str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCenterCoord/lat',xml_domnode));
    output_meta.GeoData.SCP.LLH.Lon=str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCenterCoord/lon',xml_domnode));
    output_meta.GeoData.SCP.LLH.HAE=str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneAverageHeight',xml_domnode));
end
ecf=geodetic_to_ecf([output_meta.GeoData.SCP.LLH.Lat output_meta.GeoData.SCP.LLH.Lon output_meta.GeoData.SCP.LLH.HAE]);
output_meta.GeoData.SCP.ECF.X=ecf(1);
output_meta.GeoData.SCP.ECF.Y=ecf(2);
output_meta.GeoData.SCP.ECF.Z=ecf(3);
if isempty(geo_xml) % If georeference data is available, we will just compute this later
    output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[1]/lat',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[1]/lon',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[3]/lat',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[3]/lon',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[4]/lat',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[4]/lon',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[2]/lat',xml_domnode));
    output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=...
        str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCornerCoord[2]/lon',xml_domnode));
end

%% Grid
switch char(xp.evaluate('level1Product/setup/orderInfo/projection',xml_domnode))
    case 'GROUNDRANGE'
        output_meta.Grid.ImagePlane='GROUND';
    case 'SLANTRANGE'
        output_meta.Grid.ImagePlane='SLANT';
    otherwise
        output_meta.Grid.ImagePlane='OTHER';
end
if strncmpi(char(xp.evaluate('level1Product/productSpecific/complexImageInfo/imageCoordinateType',xml_domnode)),'ZERODOPPLER',11)
    output_meta.Grid.Type='RGZERO';
end
ss_rg_s=str2double(xp.evaluate('level1Product/productInfo/imageDataInfo/imageRaster/rowSpacing',xml_domnode));
output_meta.Grid.Row.SS=ss_rg_s*SPEED_OF_LIGHT/2;
% Other ways to get values very close to this computation of Row.SS:
% output_meta.Grid.Row.SS=(1/str2double(xp.evaluate('level1Product/productSpecific/complexImageInfo/commonRSF',xml_domnode)))*SPEED_OF_LIGHT/2;
% output_meta.Grid.Row.SS=str2double(xp.evaluate('level1Product/productSpecific/complexImageInfo/projectedSpacingRange/slantRange',xml_domnode));
output_meta.Grid.Col.SS=str2double(xp.evaluate('level1Product/productSpecific/complexImageInfo/projectedSpacingAzimuth',xml_domnode));
% The following are nearly equal and should also be able to provide Col.SS
% zd_vel=str2double(xp.evaluate('level1Product/processing/geometry/zeroDopplerVelocity/velocity',xml_domnode))
% output_meta.Grid.Col.SS=zd_vel*str2double(xp.evaluate('level1Product/productInfo/imageDataInfo/imageRaster/columnSpacing',xml_domnode))
% output_meta.Grid.Col.SS=zd_vel/str2double(xp.evaluate('level1Product/productSpecific/complexImageInfo/commonPRF',xml_domnode))
% How Lockheed does it:
% output_meta.Grid.Col.SS=velocity_mag_CA_SCP * ss_zd_s * DRateSFPoly(1,1);
output_meta.Grid.Row.Sgn=-1; % Always true for TSX
output_meta.Grid.Col.Sgn=-1; % Always true for TSX
fc=str2double(xp.evaluate('level1Product/instrument/radarParameters/centerFrequency',xml_domnode));
output_meta.Grid.Row.ImpRespBW=2*str2double(xp.evaluate('level1Product/processing/processingParameter/rangeLookBandwidth',xml_domnode))/SPEED_OF_LIGHT;
dop_bw=str2double(xp.evaluate('level1Product/processing/processingParameter/azimuthLookBandwidth',xml_domnode)); % Doppler bandwidth
ss_zd_s=str2double(xp.evaluate(... % Image column spacing in zero doppler time (seconds)
    'level1Product/productInfo/imageDataInfo/imageRaster/columnSpacing',...
    xml_domnode)); % Always positive, regardless of look direction
output_meta.Grid.Col.ImpRespBW=dop_bw*ss_zd_s/output_meta.Grid.Col.SS; % Convert to azimuth spatial bandwidth (cycles per meter)
output_meta.Grid.Row.KCtr=2*fc/SPEED_OF_LIGHT;
output_meta.Grid.Col.KCtr=0;
output_meta.Grid.Row.DeltaK1=-output_meta.Grid.Row.ImpRespBW/2;
output_meta.Grid.Row.DeltaK2=-output_meta.Grid.Row.DeltaK1;
% Below describes only the case of the wrapped spectrum, which isn't always
% the case. derived_sicd_fields will populate Col.DeltaK1/2 more accurately
% later.
% output_meta.Grid.Col.DeltaK1=-(1/output_meta.Grid.Col.SS)/2;
% output_meta.Grid.Col.DeltaK2=-output_meta.Grid.Col.DeltaK1;
output_meta.Grid.Row.DeltaKCOAPoly=0;
output_meta.Grid.Row.WgtType.WindowName=...
    upper(char(xp.evaluate('level1Product/processing/processingParameter/rangeWindowID',xml_domnode)));
if strcmpi(output_meta.Grid.Row.WgtType.WindowName,'HAMMING') % The usual TSX weigting
    output_meta.Grid.Row.WgtType.Parameter.name = 'COEFFICIENT';
    output_meta.Grid.Row.WgtType.Parameter.value = char(xp.evaluate(...
        'level1Product/processing/processingParameter/rangeWindowCoefficient',...
        xml_domnode));
    a = str2double(output_meta.Grid.Row.WgtType.Parameter.value); % Generalized Hamming window parameter
    % derived_sicd_fields will populate Grid.Row.WgtFunct from WgtType later
    % Numerical solution of the analytic form of the Hamming IPR
    row_broadening_factor = 2*fzero(@(x) ...
        a*(sin(pi*x)/(pi*x)) + ((1-a)*(sin(pi*(x-1))/(pi*(x-1)))/2) + ...
        ((1-a)*(sin(pi*(x+1))/(pi*(x+1)))/2) - a/sqrt(2), .1);
    output_meta.Grid.Row.ImpRespWid = row_broadening_factor/output_meta.Grid.Row.ImpRespBW;
end
output_meta.Grid.Col.WgtType.WindowName=...
    upper(char(xp.evaluate('level1Product/processing/processingParameter/azimuthWindowID',xml_domnode)));
if strcmpi(output_meta.Grid.Col.WgtType.WindowName,'HAMMING') % The usual TSX weigting
    output_meta.Grid.Col.WgtType.Parameter.name = 'COEFFICIENT';
    output_meta.Grid.Col.WgtType.Parameter.value = char(xp.evaluate(...
        'level1Product/processing/processingParameter/azimuthWindowCoefficient',...
        xml_domnode));
    a = str2double(output_meta.Grid.Row.WgtType.Parameter.value); % Generalized Hamming window parameter
    % derived_sicd_fields will populate Grid.Row.WgtFunct from WgtType later
    col_broadening_factor = 2*fzero(@(x) ...
        a*(sin(pi*x)/(pi*x)) + ((1-a)*(sin(pi*(x-1))/(pi*(x-1)))/2) + ...
        ((1-a)*(sin(pi*(x+1))/(pi*(x+1)))/2) - a/sqrt(2), .1);
    output_meta.Grid.Col.ImpRespWid = col_broadening_factor/output_meta.Grid.Col.ImpRespBW;
end
% Grid.Col.DeltaKCOAPoly handled in "per band" section

%% Timeline
output_meta.Timeline=struct(); % Placeholder
% Most subfields added below in "per band" section

%% Position
% Compute state vectors
num_state_vectors=str2double(xp.evaluate(...
    'count(level1Product/platform/orbit/stateVec)',xml_domnode));
state_vector_T  = zeros(1,num_state_vectors);
state_vector_T_frac  = zeros(1,num_state_vectors);
state_vector_X  = zeros(1,num_state_vectors);
state_vector_Y  = zeros(1,num_state_vectors);
state_vector_Z  = zeros(1,num_state_vectors);
% state_vector_VX = zeros(1,num_state_vectors);
% state_vector_VY = zeros(1,num_state_vectors);
% state_vector_VZ = zeros(1,num_state_vectors);
for i=1:num_state_vectors
    timeStamp = char(xp.evaluate(...
        ['level1Product/platform/orbit/stateVec[' num2str(i) ']/timeUTC'],...
        xml_domnode));
    [state_vector_T(i), state_vector_T_frac(i)] = datenum_w_frac(timeStamp);
    
    state_vector_X(i) = str2double(xp.evaluate(...
        ['level1Product/platform/orbit/stateVec[' num2str(i) ']/posX'],...
        xml_domnode));
    state_vector_Y(i) = str2double(xp.evaluate(...
        ['level1Product/platform/orbit/stateVec[' num2str(i) ']/posY'],...
        xml_domnode));
    state_vector_Z(i) = str2double(xp.evaluate(...
        ['level1Product/platform/orbit/stateVec[' num2str(i) ']/posZ'],...
        xml_domnode));
%     state_vector_VX(i) = str2double(xp.evaluate(...
%         ['level1Product/platform/orbit/stateVec[' num2str(i) ']/velX'],...
%         xml_domnode));
%     state_vector_VY(i) = str2double(xp.evaluate(...
%         ['level1Product/platform/orbit/stateVec[' num2str(i) ']/velY'],...
%         xml_domnode));
%     state_vector_VZ(i) = str2double(xp.evaluate(...
%         ['level1Product/platform/orbit/stateVec[' num2str(i) ']/velZ'],...
%         xml_domnode));
end

%% RadarCollection
pol_bands=str2double(xp.evaluate('count(level1Product/productComponents/imageData/polLayer)',xml_domnode));
pols=cell(pol_bands,1);
for i=1:pol_bands
    pols{i}=char(xp.evaluate(['level1Product/productComponents/imageData[@layerIndex="' num2str(i) '"]/polLayer'],xml_domnode));
    output_meta.RadarCollection.RcvChannels.ChanParameters(i).TxRcvPolarization=[pols{i}(1) ':' pols{i}(2)];
end
tx_pol = [pols{:}]; tx_pol = unique(tx_pol(1:2:end));
if numel(tx_pol)==1
    output_meta.RadarCollection.TxPolarization = tx_pol;
else
    output_meta.RadarCollection.TxPolarization = 'SEQUENCE';
    for i = 1:numel(tx_pol)
        output_meta.RadarCollection.TxSequence(i).TxStep = i;
        output_meta.RadarCollection.TxSequence(i).TxPolarization = tx_pol(i);
    end
end

%% ImageFormation
output_meta.ImageFormation.RcvChanProc=struct('NumChanProc',1,'PRFScaleFactor',1,'ChanIndex',1);
output_meta.ImageFormation.ImageFormAlgo='RMA';
output_meta.ImageFormation.TStartProc=0;
if any(strncmpi(output_meta.CollectionInfo.RadarMode.ModeID,{'SL','HS'},2)) % SL (spotlight), HS (high-resolution spotlight)
    output_meta.ImageFormation.STBeamComp='SV';
else % SM (stripmap), SC (ScanSAR)
    output_meta.ImageFormation.STBeamComp='GLOBAL';
end
output_meta.ImageFormation.ImageBeamComp='NO';
output_meta.ImageFormation.AzAutofocus='NO';
output_meta.ImageFormation.RgAutofocus='NO';
output_meta.RMA.RMAlgoType='OMEGA_K';
output_meta.RMA.ImageType='INCA';
output_meta.RMA.INCA.FreqZero=fc;

%% SCPCOA
output_meta.SCPCOA.SideOfTrack=char(xp.evaluate('level1Product/productInfo/acquisitionInfo/lookDirection',xml_domnode));
output_meta.SCPCOA.SideOfTrack=upper(output_meta.SCPCOA.SideOfTrack(1));
if output_meta.SCPCOA.SideOfTrack=='L' % Column order flipped in COS from SICD for left-looking
    ss_zd_s=-ss_zd_s;
end
% Most subfields added below in "per band" section

%% Process fields specific to each polarimetric band
band_independent_meta=output_meta; % Values that are consistent across all bands
grouped_meta=cell(pol_bands,1);
for pol_i=1:pol_bands
    output_meta=band_independent_meta;
    
    %% Timeline
    [collectStart, collectStartFrac]=datenum_w_frac(char(xp.evaluate(...
        ['level1Product/instrument/settings[polLayer="' pols{pol_i} '"]/settingRecord/dataSegment/startTimeUTC'],...
        xml_domnode)));
    [collectEnd, collectEndFrac]=datenum_w_frac(char(xp.evaluate(...
        ['level1Product/instrument/settings[polLayer="' pols{pol_i} '"]/settingRecord/dataSegment/stopTimeUTC'],...
        xml_domnode)));
    % We loose a bit of precision when assigning the SICD CollectStart
    % field, since a MATLAB serial date number just doesn't have enough
    % bits to handle the full precision given in the TSX XML. However, all
    % relative times within the SICD metadata structure will be computed at
    % full precision.
    output_meta.Timeline.CollectStart=collectStart + (collectStartFrac/SECONDS_IN_A_DAY);
    output_meta.Timeline.CollectDuration=...
        round((collectEnd-collectStart)*SECONDS_IN_A_DAY) + ... % Convert days to seconds
        (collectEndFrac-collectStartFrac); % Handle fractional seconds
    output_meta.Timeline.IPP.Set.TStart=0;
    output_meta.Timeline.IPP.Set.TEnd=output_meta.Timeline.CollectDuration;
    output_meta.Timeline.IPP.Set.IPPStart=uint32(0);
    % Two options here to assure consistency within the IPP structure.
    % Either way, we don't understand why these two sources of metadata
    % don't exactly match up, although they are really close.
    % 1) Assume numberOfRows is the number of IPPs and derive prf/IPPPoly
    % (under the assumption the included PRF is just approximate):
    % output_meta.Timeline.IPP.Set.IPPEnd=uint32(str2double(char(xp.evaluate(...
    %     ['level1Product/instrument/settings[polLayer="' pols{pol_i} ...
    %     '"]/settingRecord/dataSegment/numberOfRows'],xml_domnode))) ...
    %     - 1); % Since these indices are zero-based.
    % output_meta.Timeline.IPP.Set.IPPPoly = [0; ...
    %     double(output_meta.Timeline.IPP.Set.IPPEnd)/output_meta.Timeline.IPP.Set.TEnd];
    % 2) Assume prf is exact and derive number of IPPs (because we don't
    % really understand exactly what numberOfRows is):
    prf=str2double(xp.evaluate(...
        ['level1Product/instrument/settings[polLayer="' pols{pol_i} '"]/settingRecord/PRF'],...
        xml_domnode));
    output_meta.Timeline.IPP.Set.IPPPoly=[0; prf];
    output_meta.Timeline.IPP.Set.IPPEnd = ...
        uint32(polyval(output_meta.Timeline.IPP.Set.IPPPoly(end:-1:1),...
        double(output_meta.Timeline.IPP.Set.TEnd)));
    
    %% Position
    % Polynomial is computed with respect to time from start of collect
    state_vector_T_band = round((state_vector_T-collectStart)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (state_vector_T_frac-collectStartFrac); % Handle fractional seconds
    % sv2poly.m shows ways to determine best polynomial order, but 5th is almost always best
    polyorder=min(5, numel(state_vector_T_band) - 1);
    P_x  = polyfit(state_vector_T_band, state_vector_X, polyorder);
    P_y  = polyfit(state_vector_T_band, state_vector_Y, polyorder);
    P_z  = polyfit(state_vector_T_band, state_vector_Z, polyorder);
    output_meta.Position.ARPPoly.X  = P_x(end:-1:1).';
    output_meta.Position.ARPPoly.Y  = P_y(end:-1:1).';
    output_meta.Position.ARPPoly.Z  = P_z(end:-1:1).';

    %% RadarCollection
    % output_meta.RadarCollection.RefFreqIndex=uint32(0); % Absence of this field means all frequencies are true values
    bw=str2double(xp.evaluate(['level1Product/instrument/settings[polLayer="' pols{pol_i} '"]/rxBandwidth'],xml_domnode));
    output_meta.RadarCollection.TxFrequency.Min=fc-(bw/2); % fc calculated in grid section
    output_meta.RadarCollection.TxFrequency.Max=fc+(bw/2);
    % Pulse length scaling factor from personal communication with Juergen Janoth, Head of Application Development, Infoterra
    output_meta.RadarCollection.Waveform.WFParameters.TxPulseLength=str2double(xp.evaluate(...
        'level1Product/processing/processingParameter/rangeCompression/chirps/referenceChirp/pulseLength',...
        xml_domnode)) * 32 / (3.29658384e8);
    output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth=bw; % This is actually receive BW.  We assume its the same.
    output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart=...
        output_meta.RadarCollection.TxFrequency.Min;
    output_meta.RadarCollection.Waveform.WFParameters.TxFMRate=...
        bw/output_meta.RadarCollection.Waveform.WFParameters.TxPulseLength;
    output_meta.RadarCollection.Waveform.WFParameters.RcvDemodType='CHIRP';
    output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate=...
        str2double(xp.evaluate(['level1Product/instrument/settings[polLayer="' pols{pol_i} '"]/RSF'],xml_domnode));
    output_meta.RadarCollection.Waveform.WFParameters.RcvWindowLength=...
        str2double(xp.evaluate(['level1Product/instrument/settings[polLayer="' pols{pol_i} '"]/settingRecord/echowindowLength'],xml_domnode))/...
        output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate;
    output_meta.RadarCollection.Waveform.WFParameters.RcvFMRate=0; % True for RcvDemodType='CHIRP'
    
    %% ImageFormation
    output_meta.ImageFormation.RcvChanProc.ChanIndex=pol_i;
    output_meta.ImageFormation.TEndProc=output_meta.Timeline.CollectDuration;
    output_meta.ImageFormation.TxFrequencyProc.MinProc=...
        output_meta.RadarCollection.TxFrequency.Min;
    output_meta.ImageFormation.TxFrequencyProc.MaxProc=...
        output_meta.RadarCollection.TxFrequency.Max;
    output_meta.ImageFormation.TxRcvPolarizationProc=...
        output_meta.RadarCollection.RcvChannels.ChanParameters(pol_i).TxRcvPolarization;

    %% RMA
    if ~isempty(geo_xml) % Use geoReference file which is much more precise
        [t_ref, t_ref_frac] = datenum_w_frac(char(xp.evaluate(...
            'geoReference/geolocationGrid/gridReferenceTime/tReferenceTimeUTC', geo_xml)));
        az_offset = round((t_ref-collectStart)*SECONDS_IN_A_DAY) + ... % Convert to seconds
            (t_ref_frac-collectStartFrac); % Handle fractional seconds
        TimeCA_scp = str2double(char(xp.evaluate(... % Zero doppler time for SCP
            [scp_geo_node_str '/t'], geo_xml)));
        az_shift = 0;
        for j = 1:str2double(xp.evaluate('count(geoReference/signalPropagationEffects/azimuthShift)',geo_xml))
            az_shift = az_shift + ...
            str2double(xp.evaluate(['geoReference/signalPropagationEffects/azimuthShift[' num2str(j) ']/coefficient'], geo_xml));
        end
        output_meta.RMA.INCA.TimeCAPoly = TimeCA_scp + az_offset - az_shift; % Relative to start of collect, not geogrid reference time
        azimuth_time_scp = t_ref;
        azimuth_time_frac = t_ref_frac + TimeCA_scp; % Might actually be >1, but that's OK
        range_time_scp = ...
            str2double(xp.evaluate('geoReference/geolocationGrid/gridReferenceTime/tauReferenceTime', geo_xml)) + ...
            str2double(xp.evaluate([scp_geo_node_str '/tau'], geo_xml));
        rg_delay = 0;
        for j = 1:str2double(xp.evaluate('count(geoReference/signalPropagationEffects/rangeDelay)',geo_xml))
            rg_delay = rg_delay + ...
            str2double(xp.evaluate(['geoReference/signalPropagationEffects/rangeDelay[' num2str(j) ']/coefficient'], geo_xml));
        end
        output_meta.RMA.INCA.R_CA_SCP = (range_time_scp-rg_delay)*SPEED_OF_LIGHT/2;
    else % Main TSX XML also defines its own center point
        [azimuth_time_scp, azimuth_time_frac] = datenum_w_frac(char(xp.evaluate(... % Zero doppler time for SCP
            'level1Product/productInfo/sceneInfo/sceneCenterCoord/azimuthTimeUTC',...
            xml_domnode)));
        output_meta.RMA.INCA.TimeCAPoly = ...
            round((azimuth_time_scp-collectStart)*SECONDS_IN_A_DAY) + ... % Convert to seconds
            (azimuth_time_frac-collectStartFrac); % Handle fractional seconds
        range_time_scp=str2double(xp.evaluate(...
            'level1Product/productInfo/sceneInfo/sceneCenterCoord/rangeTime',...
            xml_domnode)); % Range for SCP
        output_meta.RMA.INCA.R_CA_SCP=range_time_scp*SPEED_OF_LIGHT/2; % Range for SCP
    end
    output_meta.RMA.INCA.TimeCAPoly(2,1)=ss_zd_s/output_meta.Grid.Col.SS; % Convert zero doppler spacing from sec/pixels to sec/meters
    [output_meta.RMA.INCA.DopCentroidPoly, ...
        output_meta.Grid.TimeCOAPoly] = ...
        computeDopCentroidPolys(xml_domnode,pols{pol_i},...
        azimuth_time_scp,azimuth_time_frac,... % SCP zero doppler time
        range_time_scp); % SCP range time
    if strncmpi(output_meta.CollectionInfo.RadarMode.ModeID,'ST',2)
        % DopCentroidPoly only used for Stripmap and Dynamic Stripmap
        % collections.
        output_meta.RMA.INCA = rmfield(output_meta.RMA.INCA,'DopCentroidPoly');
        output_meta.Grid.Col.DeltaKCOAPoly=0;
        % This isn't exactly right.  TimeCOAPoly should not depend on which
        % SCP pixel was chosen:
        output_meta.Grid.TimeCOAPoly=output_meta.RMA.INCA.TimeCAPoly(1,1);
    else
        output_meta.RMA.INCA.DopCentroidCOA=true;
        output_meta.Grid.Col.DeltaKCOAPoly=...
            output_meta.RMA.INCA.DopCentroidPoly*ss_zd_s/output_meta.Grid.Col.SS;
    end
    
    % Compute DRateSFPoly
    % For the purposes of the DRateSFPoly computation, we ignore any
    % changes in velocity or doppler rate over the azimuth dimension.
    pos_coefs = [P_x(:) P_y(:) P_z(:)];
    % Velocity is derivate of position.
    vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
    vel_x = polyval(vel_coefs(:,1), output_meta.RMA.INCA.TimeCAPoly(1));
    vel_y = polyval(vel_coefs(:,2), output_meta.RMA.INCA.TimeCAPoly(1));
    vel_z = polyval(vel_coefs(:,3), output_meta.RMA.INCA.TimeCAPoly(1));
    vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
    r_ca = [output_meta.RMA.INCA.R_CA_SCP; 1]; % Polynomial representing range as a function of range distance from SCP
    % For the DRateSFPoly computations here (and the
    % computeDopCentroidPolys subfunction below), we use a single doppler
    % rate polynomial.  For most collects, only 3 doppler rate polynomials
    % are provided and they don't vary significantly, so this is probably
    % OK. For staring spotlight mode collects though, many more polynomials
    % are given and using a single polynomial to represent the entire
    % collect may not be sufficient.
    num_dop_rates=str2double(xp.evaluate(...
        'count(level1Product/processing/geometry/dopplerRate)',...
        xml_domnode));
    center_dop_rate=num2str(ceil(num_dop_rates/2)); % Pick the center
    num_dop_rate_coefs=str2double(xp.evaluate(...
        ['count(level1Product/processing/geometry/dopplerRate[' ...
        center_dop_rate ']/dopplerRatePolynomial/coefficient)'],...
        xml_domnode));
    dop_rate_poly = zeros(num_dop_rate_coefs,1);
    for i=1:num_dop_rate_coefs
        dop_rate_poly(i)=str2double(xp.evaluate(... % Doppler FM rate with regard to raw time
            ['level1Product/processing/geometry/dopplerRate[' ...
            center_dop_rate ']/dopplerRatePolynomial/coefficient[@exponent="' num2str(i-1) '"]'],...
            xml_domnode));
    end
    dop_rate_ref = str2double(xp.evaluate(... % Doppler FM rate reference time
        ['level1Product/processing/geometry/dopplerRate[' ...
        center_dop_rate ']/dopplerRatePolynomial/referencePoint'],...
        xml_domnode));
    % Shift 1D polynomial to account for SCP
    dop_rate_poly_rg_shifted=polyshift(dop_rate_poly, range_time_scp-dop_rate_ref);
    % Scale 1D polynomial to from Hz/s^n to Hz/m^n
    dop_rate_poly_rg_scaled=dop_rate_poly_rg_shifted.*...
        (2/SPEED_OF_LIGHT).^(0:(length(dop_rate_poly)-1)).';
    output_meta.RMA.INCA.DRateSFPoly = - conv(dop_rate_poly_rg_scaled,r_ca) * ... % Multiplication of two polynomials is just a convolution of their coefficients
        SPEED_OF_LIGHT / (2 * fc * vm_ca_sq(1)); % Assumes a SGN of -1

    %% GeoData
    % Now that sensor model fields have been populated, we can populate
    % GeoData.SCP more precisely.
    ecf = point_image_to_ground([output_meta.ImageData.SCPPixel.Row;output_meta.ImageData.SCPPixel.Col],output_meta);
    output_meta.GeoData.SCP.ECF.X=ecf(1);
    output_meta.GeoData.SCP.ECF.Y=ecf(2);
    output_meta.GeoData.SCP.ECF.Z=ecf(3);
    llh=ecf_to_geodetic([output_meta.GeoData.SCP.ECF.X output_meta.GeoData.SCP.ECF.Y output_meta.GeoData.SCP.ECF.Z]);
    output_meta.GeoData.SCP.LLH.Lat=llh(1);
    output_meta.GeoData.SCP.LLH.Lon=llh(2);
    output_meta.GeoData.SCP.LLH.HAE=llh(3);
    
    %% SCPCOA
    output_meta = derived_sicd_fields(output_meta);
    % Can also use incidence angle from the TSX XML metdata
    % output_meta.SCPCOA.IncidenceAng=str2double(xp.evaluate('level1Product/productInfo/sceneInfo/sceneCenterCoord/incidenceAngle',xml_domnode));
    % output_meta.SCPCOA.GrazeAng=90-output_meta.SCPCOA.IncidenceAng;
    % Could use XML velocity values directly
    % P_vx = polyfit(state_vector_T_band, state_vector_VX, polyorder);
    % P_vy = polyfit(state_vector_T_band, state_vector_VY, polyorder);
    % P_vz = polyfit(state_vector_T_band, state_vector_VZ, polyorder);
    % output_meta.SCPCOA.ARPVel.X = polyval(P_vx,output_meta.SCPCOA.SCPTime);
    % output_meta.SCPCOA.ARPVel.Y = polyval(P_vy,output_meta.SCPCOA.SCPTime);
    % output_meta.SCPCOA.ARPVel.Z = polyval(P_vz,output_meta.SCPCOA.SCPTime);

    grouped_meta{pol_i}=output_meta;
end
output_meta=grouped_meta; % Cell array with metadata struct for each band

end


function [DopCentroidPoly, TimeCOAPoly] = ...
        computeDopCentroidPolys(domnode, polarization, ...
        t_utc_scp, t_utc_scp_frac, t_rg_scp) % SCP info with which these polynomials will be with to respect to
% From the paper "TerraSAR-X Deskew Description", Michael Stewart,
% December 11, 2008

SECONDS_IN_A_DAY = 24*60*60;

xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

% Read input parameters from metadata
N=str2double(xp.evaluate(...
    ['count(level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate)'],...
    domnode));
M=49; % Tunable parameter.  Set to 49 in the paper.
t_utc_ref=char(xp.evaluate(... % Reference UTC time for doppler centroid polynomial fit
    ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(ceil(N/2)) ']/timeUTC'],...
    domnode));
[t_utc_ref, t_utc_ref_frac] = datenum_w_frac(t_utc_ref);
t_rg_ref=str2double(xp.evaluate(... % Reference range time for doppler centroid polynomial fit
    ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(ceil(N/2)) ']/combinedDoppler/referencePoint'],...
    domnode));
[t_utc_start,t_utc_start_frac]=datenum_w_frac(char(xp.evaluate(...
        ['level1Product/instrument/settings[polLayer="' upper(polarization) '"]/settingRecord/dataSegment/startTimeUTC'],...
        domnode)));
[det_raw,det_raw_frac,det_rg_max,det_rg_min,de_c0,de_c1]=deal(zeros(N,1));
for i=1:N
    det_raw_str=char(xp.evaluate(... % Raw times of doppler estimates.  We use this for center of aperture time.
        ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(i) ']/timeUTC'],...
        domnode));
    [det_raw(i), det_raw_frac(i)]=datenum_w_frac(det_raw_str);
    det_rg_max(i)=str2double(xp.evaluate(... % Validity range max of doppler estimates
        ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(i) ']/combinedDoppler/validityRangeMax'],...
        domnode));
    det_rg_min(i)=str2double(xp.evaluate(... % Validity range min of doppler estimates
        ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(i) ']/combinedDoppler/validityRangeMin'],...
        domnode));
    de_c0(i)=str2double(xp.evaluate(... % Constant terms of the doppler centroid estimate of valid range times
        ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(i) ']/combinedDoppler/coefficient[@exponent="0"]'],...
        domnode));
    de_c1(i)=str2double(xp.evaluate(... % Linear terms of the doppler centroid estimate of valid range times
        ['level1Product/processing/doppler/dopplerCentroid[polLayer="' upper(polarization) '"]/dopplerEstimate[' num2str(i) ']/combinedDoppler/coefficient[@exponent="1"]'],...
        domnode));
end
num_dop_rates=str2double(xp.evaluate(...
    'count(level1Product/processing/geometry/dopplerRate)',...
    domnode));
center_dop_rate=num2str(ceil(num_dop_rates/2)); % Pick the center
fm_dop=str2double(xp.evaluate(... % Doppler FM rate with regard to raw time
    ['level1Product/processing/geometry/dopplerRate[' center_dop_rate ']/dopplerRatePolynomial/coefficient[@exponent="0"]'],...
    domnode)); % Varies only slightly across the image.  We use a single value here.
ss_zd_s=str2double(xp.evaluate(... % Image column spacing in zero doppler time
    'level1Product/productInfo/imageDataInfo/imageRaster/columnSpacing',...
    domnode));
lookdir=char(xp.evaluate('level1Product/productInfo/acquisitionInfo/lookDirection',domnode));
ss_zd_m=str2double(xp.evaluate('level1Product/productSpecific/complexImageInfo/projectedSpacingAzimuth',domnode));

% Doppler Centroid Grid Evaluation (section 2.2 in paper)
% Compute differential range and zero doppler times for each of the TSX
% discrete doppler centroid estimates
[diff_t_raw,diff_t_rg,f_dc,diff_t_zd]=deal(zeros(M,N));
for n=1:N
    delta_t_rg=(det_rg_max(n)-det_rg_min(n))/(M-1);
    for m=1:M
        diff_t_raw(m,n)=round((det_raw(n)-t_utc_ref)*SECONDS_IN_A_DAY) + ... % Convert from days to seconds
            (det_raw_frac(n)-t_utc_ref_frac);
        diff_t_rg(m,n)=det_rg_min(n)-t_rg_ref+((m-1)*delta_t_rg);
        f_dc(m,n)=de_c0(n)+de_c1(n)*diff_t_rg(m,n); % Sampled values of Doppler centroid
        diff_t_zd(m,n)=diff_t_raw(m,n)-(f_dc(m,n)/fm_dop); % TimeCA relative to reference time
    end
end
% diff_t_raw is center-of-aperture time relative to reference time.  We
% want it from start of collect.
t_coa = diff_t_raw + round((t_utc_ref - t_utc_start)*SECONDS_IN_A_DAY) + ...
    (t_utc_ref_frac - t_utc_start_frac);

% The paper computes a polynomial as a function of differential time.  For
% SICD, we need polynomial as a function of range/azimuth distance from SCP
% in meters.
range_scp_m = (diff_t_rg + t_rg_ref - t_rg_scp) * (SPEED_OF_LIGHT/2);
azimuth_scp_m = ss_zd_m * (diff_t_zd + ...
    round((t_utc_ref - t_utc_scp) * SECONDS_IN_A_DAY) + ...
    (t_utc_ref_frac - t_utc_scp_frac)) / ss_zd_s;
if lookdir(1)=='L' % SICD is in view-from-above not oriented by collection direction
    azimuth_scp_m=-azimuth_scp_m;
end

% Least squares fit (section 2.3 in paper)
% Compute the following continuous function for doppler centroid
% f_dc(azimuth_scp_m,range_scp_m)=...
%    x(1,1) + x(2,1)*range_scp_m + x(3,1)*range_scp_m*range_scp_m+...
%    x(1,2)*azimuth_scp_m + x(1,3)*azimuth_scp_m*azimuth_scp_m+...
%    x(2,2)*azimuth_scp_m*range_scp_m;
% A*x = b
% We add x(2,3),x(3,2), and x(3,3), even though they were not in the paper.
a=[ones(M*N,1) ...
   range_scp_m(:) ...
   range_scp_m(:).^2 ...
   azimuth_scp_m(:) ...
   range_scp_m(:).*azimuth_scp_m(:) ...
   (range_scp_m(:).^2).*azimuth_scp_m(:) ...
   azimuth_scp_m(:).^2 ...
   range_scp_m(:).*(azimuth_scp_m(:).^2) ...
   (range_scp_m(:).^2).*(azimuth_scp_m(:).^2)];
[b_dc,b_coa]=deal(zeros(size(a,2),1));
for i=1:size(a,2)
   b_dc(i)=sum(f_dc(:).*a(:,i)); % Doppler centroid
   b_coa(i)=sum(t_coa(:).*a(:,i)); % Center of aperture time
end
A=zeros(size(a,2));
for i=1:size(a,2)
    for j=1:size(a,2)
        A(i,j)=sum(a(:,i).*a(:,j));
    end
end
old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
x=A\b_dc; % MATLAB often flags this as badly scaled, but results still appear valid
x2=A\b_coa;
warning(old_warning_state);
DopCentroidPoly=reshape(x,3,3); % Polynomial coefficients for computing the doppler
                                % centroid as a continous function of range
                                % and azimuth distance from SCP.
TimeCOAPoly=reshape(x2,3,3);    % Polynomial coefficients for computing
                                % coa time as a continous function of
                                % range and azimuth distance from SCP.

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
    if isnan(datenum_frac), datenum_frac = 0; end;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////