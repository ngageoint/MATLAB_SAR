function [ output_meta ] = meta2sicd_rs_xml( xml_domnode, beta_domnode, noise_domnode )
%META2SICD_RS_XML Converts Radarsat product.xml description into a SICD-style metadata structure
%
% Takes as input a Document Object Model (DOM) node from the RS2
% product.xml descriptor file.
%
% This function does NOT handle ScanSAR datasets.
%
% There is an outstanding question with regard to the meaning of the
% pulseRepetitionFrequency as provided.  See comments in Timeline section
% below.
%
% Written by: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup
SECONDS_IN_A_DAY = 24*60*60;
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

%% CollectionInfo
output_meta.CollectionInfo.CollectorName=char(xp.evaluate(...
    xpath_str({'product','sourceAttributes','satellite'}),xml_domnode));
if strcmpi(output_meta.CollectionInfo.CollectorName, 'RADARSAT-2')
    gen = 'RS2';
elseif strncmpi(output_meta.CollectionInfo.CollectorName, 'RCM', 3)
    gen = 'RCM';
end
[rawDataStartTime, rawDataStartTimeFrac] = datenum_w_frac(char(xp.evaluate(...
    xpath_str({'product','sourceAttributes','rawDataStartTime'}),...
    xml_domnode)));
if strcmp(gen,'RS2')
    output_meta.CollectionInfo.CoreName = [... % Start with NGA-like prefix
        upper(datestr(rawDataStartTime,'ddmmmyy')) 'RS02' ...
        char(xp.evaluate(xpath_str({'product','sourceAttributes','imageId'}),xml_domnode))];
elseif strcmp(gen, 'RCM')
    output_meta.CollectionInfo.CoreName = [... % Start with NGA-like prefix
        upper(datestr(rawDataStartTime,'ddmmmyy')) ...
        'RCM' output_meta.CollectionInfo.CollectorName(end) ...
        datestr(rawDataStartTime,'HHMMSS')]; % Make time of day unique identier within day
    % ScanSAR might need multiple CoreNames, one for each burst?
end
output_meta.CollectionInfo.CollectType='MONOSTATIC';
output_meta.CollectionInfo.RadarMode.ModeID=char(xp.evaluate(...
    xpath_str({'product','sourceAttributes','beamModeMnemonic'}),xml_domnode));
beamMode = char(xp.evaluate(xpath_str({...
    'product','sourceAttributes','beamMode'}),xml_domnode));
acqType = char(xp.evaluate(xpath_str({...
    'product','sourceAttributes','radarParameters','acquisitionType'}),xml_domnode));
if ((~isempty(beamMode) && strncmpi(beamMode, 'SPOTLIGHT', 9)) || ...
   (~isempty(acqType) && strncmpi(acqType, 'SPOTLIGHT', 9)) || ...
   ~isempty(strfind(output_meta.CollectionInfo.RadarMode.ModeID,'SL')))
    output_meta.CollectionInfo.RadarMode.ModeType = 'SPOTLIGHT';
elseif strcmpi(output_meta.CollectionInfo.RadarMode.ModeID(1:2), 'SC')
    error('META2SICD_RS_XML:RS_SCANSAR', 'ScanSAR mode data is not currently handled.');
else % Finally assume it's stripmap
    output_meta.CollectionInfo.RadarMode.ModeType = 'STRIPMAP';
end
if strcmp(gen, 'RS2') % All RS2 data is unclassified
    output_meta.CollectionInfo.Classification='UNCLASSIFIED';
elseif strcmp(gen, 'RCM') % RCM has this as a specific field
    output_meta.CollectionInfo.Classification=upper(char(xp.evaluate(...
        xpath_str({'product', 'securityAttributes', 'securityClassification'}), xml_domnode)));
end

%% ImageCreation
output_meta.ImageCreation.Application=char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','generalProcessingInformation','softwareVersion'}),...
    xml_domnode));
output_meta.ImageCreation.DateTime=datenum(char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','generalProcessingInformation','processingTime'}),...
    xml_domnode)),'yyyy-mm-ddTHH:MM:SS.FFF');
output_meta.ImageCreation.Site=char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','generalProcessingInformation','processingFacility'}),...
    xml_domnode));
output_meta.ImageCreation.Profile='Prototype';

%% ImageData
if strcmp(gen, 'RS2')
    output_meta.ImageData.NumRows=uint32(str2double(xp.evaluate(...
       xpath_str({'product','imageAttributes','rasterAttributes','numberOfSamplesPerLine'}),...
       xml_domnode)));
    output_meta.ImageData.NumCols=uint32(str2double(xp.evaluate(...
       xpath_str({'product','imageAttributes','rasterAttributes','numberOfLines'}),...
        xml_domnode)));
elseif strcmp(gen, 'RCM')
    output_meta.ImageData.NumRows=uint32(str2double(xp.evaluate(...
       xpath_str({'product','sceneAttributes','imageAttributes','samplesPerLine'}),...
       xml_domnode)));
    output_meta.ImageData.NumCols=uint32(str2double(xp.evaluate(...
       xpath_str({'product','sceneAttributes','imageAttributes','numLines'}),...
        xml_domnode)));
end
output_meta.ImageData.FullImage=output_meta.ImageData;
output_meta.ImageData.FirstRow=uint32(0); output_meta.ImageData.FirstCol=uint32(0);
output_meta.ImageData.PixelType='RE16I_IM16I';  % RS2 always 16-bit
if strcmp(gen, 'RCM') && str2double(xp.evaluate(...  % RCM can be 16 or 32
       xpath_str({'product','imageReferenceAttributes','rasterAttributes','bitsPerSample'}),...
        xml_domnode)) == 32
    output_meta.ImageData.PixelType='RE32F_IM32F';
end
% Seems that all pixels are always valid
output_meta.ImageData.ValidData.Vertex(1).Row = uint32(0);
output_meta.ImageData.ValidData.Vertex(1).Col = uint32(0);
output_meta.ImageData.ValidData.Vertex(2).Row = uint32(0);
output_meta.ImageData.ValidData.Vertex(2).Col = output_meta.ImageData.NumCols-1;
output_meta.ImageData.ValidData.Vertex(3).Row = output_meta.ImageData.NumRows-1;
output_meta.ImageData.ValidData.Vertex(3).Col = output_meta.ImageData.NumCols-1;
output_meta.ImageData.ValidData.Vertex(4).Row = output_meta.ImageData.NumRows-1;
output_meta.ImageData.ValidData.Vertex(4).Col = uint32(0);

%% SCP
if strcmp(gen, 'RS2')
    im_at_str = 'imageAttributes';
elseif strcmp(gen, 'RCM')
    im_at_str = 'imageReferenceAttributes';
end
% There are many different equally valid options for picking the SCP point.
% One way is to chose the tie point that is closest to the image center.
num_tie_points=str2double(xp.evaluate(...
    ['count(' xpath_str({'product',im_at_str,'geographicInformation','geolocationGrid','imageTiePoint'}) ')'],...
    xml_domnode));
tiePointPixels = zeros(2,num_tie_points);
tiePointGeo    = zeros(3,num_tie_points);
for i=1:num_tie_points
    tiePointPixels(1,i) = str2double(xp.evaluate(...
            [xpath_str({'product',im_at_str,'geographicInformation','geolocationGrid','imageTiePoint'})...
            '[' num2str(i) ']' xpath_str({'imageCoordinate','pixel'})],...
            xml_domnode));
    tiePointPixels(2,i) = str2double(xp.evaluate(...
            [xpath_str({'product',im_at_str,'geographicInformation','geolocationGrid','imageTiePoint'})...
            '[' num2str(i) ']' xpath_str({'imageCoordinate','line'})],...
            xml_domnode));
    tiePointGeo(1,i) = str2double(xp.evaluate(...
            [xpath_str({'product',im_at_str,'geographicInformation','geolocationGrid','imageTiePoint'})...
            '[' num2str(i) ']' xpath_str({'geodeticCoordinate','latitude'})],...
            xml_domnode));
    tiePointGeo(2,i) = str2double(xp.evaluate(...
            [xpath_str({'product',im_at_str,'geographicInformation','geolocationGrid','imageTiePoint'})...
            '[' num2str(i) ']' xpath_str({'geodeticCoordinate','longitude'})],...
            xml_domnode));
    tiePointGeo(3,i) = str2double(xp.evaluate(...
            [xpath_str({'product',im_at_str,'geographicInformation','geolocationGrid','imageTiePoint'})...
            '[' num2str(i) ']' xpath_str({'geodeticCoordinate','height'})],...
            xml_domnode));
end
% Pick tie point closest to center for SCP
centerPoint = [double(output_meta.ImageData.NumRows-1)/2.0; ...
               double(output_meta.ImageData.NumCols-1)/2.0];
D = (tiePointPixels - repmat(centerPoint,1,num_tie_points));
[C,scp_index] = min( sqrt( sum(D.^2) ) );
output_meta.ImageData.SCPPixel.Row = uint32(tiePointPixels(1,scp_index));
output_meta.ImageData.SCPPixel.Col = uint32(tiePointPixels(2,scp_index));

% Sometimes lines up with SCP in Lockheed SICDs (and sometimes not):
% output_meta.ImageData.SCPPixel.Col = ceil(output_meta.ImageData.NumCols/2)-1;
% output_meta.ImageData.SCPPixel.Row = floor(output_meta.ImageData.NumRows/2);

%% GeoData
% All RS2 and RCM data we know use the WGS84 model, although it is stated
% in slightly different XML fields in the RS2 and RCM product.xml
output_meta.GeoData.EarthModel='WGS_84';
% Initially, we just seed this with a rough value.  Later we will put in
% something more precise.
output_meta.GeoData.SCP.LLH.Lat = tiePointGeo(1,scp_index);
output_meta.GeoData.SCP.LLH.Lon = tiePointGeo(2,scp_index);
output_meta.GeoData.SCP.LLH.HAE = tiePointGeo(3,scp_index);
pos_ecf = geodetic_to_ecf(tiePointGeo(:,scp_index));
output_meta.GeoData.SCP.ECF.X = pos_ecf(1);
output_meta.GeoData.SCP.ECF.Y = pos_ecf(2);
output_meta.GeoData.SCP.ECF.Z = pos_ecf(3);
% Corner coordinates will be computed later by derived_sicd_fields.
% Corners
% min_row=min(tiePointPixels(1,:));
% min_col=min(tiePointPixels(2,:));
% max_row=max(tiePointPixels(1,:));
% max_col=max(tiePointPixels(2,:));
% ll_index=find((tiePointPixels(1,:)==min_row)&(tiePointPixels(2,:)==min_col), 1);
% ul_index=find((tiePointPixels(1,:)==max_row)&(tiePointPixels(2,:)==min_col), 1);
% ur_index=find((tiePointPixels(1,:)==max_row)&(tiePointPixels(2,:)==max_col), 1);
% lr_index=find((tiePointPixels(1,:)==min_row)&(tiePointPixels(2,:)==max_col), 1);
% output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=tiePointGeo(1,ul_index);
% output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=tiePointGeo(2,ul_index);
% output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=tiePointGeo(1,ur_index);
% output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=tiePointGeo(2,ur_index);
% output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=tiePointGeo(1,lr_index);
% output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=tiePointGeo(2,lr_index);
% output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=tiePointGeo(1,ll_index);
% output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=tiePointGeo(2,ll_index);

%% State vectors
num_state_vectors=str2double(xp.evaluate(...
    ['count(' xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'}) ')'],...
    xml_domnode));
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
        [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
        '[' num2str(i) ']' xpath_str({'timeStamp'})],...
        xml_domnode));
    [state_vector_T(i), state_vector_T_frac(i)] = datenum_w_frac(timeStamp);
    
    state_vector_X(i) = str2double(xp.evaluate(...
        [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
        '[' num2str(i) ']' xpath_str({'xPosition'})],...
        xml_domnode));
    state_vector_Y(i) = str2double(xp.evaluate(...
        [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
        '[' num2str(i) ']' xpath_str({'yPosition'})],...
        xml_domnode));
    state_vector_Z(i) = str2double(xp.evaluate(...
        [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
        '[' num2str(i) ']' xpath_str({'zPosition'})],...
        xml_domnode));
%     state_vector_VX(i) = str2double(xp.evaluate(...
%         [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
%         '[' num2str(i) ']' xpath_str({'xVelocity'})],...
%         xml_domnode));
%     state_vector_VY(i) = str2double(xp.evaluate(...
%         [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
%         '[' num2str(i) ']' xpath_str({'yVelocity'})],...
%         xml_domnode));
%     state_vector_VZ(i) = str2double(xp.evaluate(...
%         [xpath_str({'product','sourceAttributes','orbitAndAttitude','orbitInformation','stateVector'})...
%         '[' num2str(i) ']' xpath_str({'zVelocity'})],...
%         xml_domnode));
end
state_vector_T = round((state_vector_T-rawDataStartTime)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (state_vector_T_frac-rawDataStartTimeFrac); % Handle fractional seconds
% sv2poly.m shows ways to determine best polynomial order, but 5th is almost always best
polyorder = min(5, numel(state_vector_T) - 1);
P_x = polyfit(state_vector_T, state_vector_X, polyorder);
P_y = polyfit(state_vector_T, state_vector_Y, polyorder);
P_z = polyfit(state_vector_T, state_vector_Z, polyorder);
output_meta.Position.ARPPoly.X = P_x(end:-1:1).';
output_meta.Position.ARPPoly.Y = P_y(end:-1:1).';
output_meta.Position.ARPPoly.Z = P_z(end:-1:1).';

%% Grid
if strcmp(char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','generalProcessingInformation','productType'}),...
    xml_domnode)),'SLC')
    output_meta.Grid.ImagePlane = 'SLANT';
    output_meta.Grid.Type = 'RGZERO';
else
    output_meta.Grid.ImagePlane = 'GROUND';
end
output_meta.Grid.Row.SS = str2double(xp.evaluate(...
    xpath_str({'product',im_at_str,'rasterAttributes','sampledPixelSpacing'}),...
    xml_domnode));
% Col.SS is derived after DRateSFPoly below, rather than used from this
% given field, so that SICD metadata can be internally consistent.
% output_meta.Grid.Col.SS = str2double(xp.evaluate(...
%     xpath_str({'product',im_at_str,'rasterAttributes','sampledLineSpacing'}),...
%     xml_domnode));
output_meta.Grid.Row.Sgn = -1; % Always true for RS2
output_meta.Grid.Col.Sgn = -1; % Always true for RS2
fc = str2double(xp.evaluate(...
    xpath_str({'product','sourceAttributes','radarParameters','radarCenterFrequency'}),...
    xml_domnode)); % Center frequency
output_meta.Grid.Row.ImpRespBW = 2*str2double(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','sarProcessingInformation','totalProcessedRangeBandwidth'}),...
    xml_domnode))/SPEED_OF_LIGHT;
dop_bw = str2double(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','sarProcessingInformation','totalProcessedAzimuthBandwidth'}),...
    xml_domnode)); % Doppler bandwidth
[zd_last, zd_last_frac] = datenum_w_frac(char(xp.evaluate(... 
    xpath_str({'product','imageGenerationParameters','sarProcessingInformation','zeroDopplerTimeLastLine'}),...
    xml_domnode)));
[zd_first, zd_first_frac] = datenum_w_frac(char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','sarProcessingInformation','zeroDopplerTimeFirstLine'}),...
    xml_domnode)));
ss_zd_s = abs(round((zd_last-zd_first)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (zd_last_frac-zd_first_frac))/... % Handle fractional seconds
        double(output_meta.ImageData.NumCols-1); % Image column spacing in zero doppler time (seconds)
output_meta.Grid.Row.KCtr = 2*fc/SPEED_OF_LIGHT;
output_meta.Grid.Col.KCtr = 0;
output_meta.Grid.Row.DeltaK1 = -output_meta.Grid.Row.ImpRespBW/2;
output_meta.Grid.Row.DeltaK2 = -output_meta.Grid.Row.DeltaK1;
output_meta.Grid.Row.DeltaKCOAPoly = 0;
% Constants used to compute weighting parameters
NUM_SAMPLES = 512;
OVERSAMPLE = 1024;
output_meta.Grid.Row.WgtType.WindowName = upper(char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','sarProcessingInformation','rangeWindow','windowName'}),...
    xml_domnode)));
if strcmpi(output_meta.Grid.Row.WgtType.WindowName,'KAISER') % The usual RS2 weigting
    output_meta.Grid.Row.WgtType.Parameter.name = 'BETA';
    output_meta.Grid.Row.WgtType.Parameter.value = char(xp.evaluate(...
        xpath_str({'product','imageGenerationParameters','sarProcessingInformation','rangeWindow','windowCoefficient'}),...
        xml_domnode));
    beta_row = str2double(output_meta.Grid.Row.WgtType.Parameter.value);
    % We don't use the Mathworks Kaiser function, so we won't be dependent on the Signal Processing Toolbox
    output_meta.Grid.Row.WgtFunct = kaiser_nosptb(NUM_SAMPLES,beta_row);
    imp_resp = abs(fft(output_meta.Grid.Row.WgtFunct, round(NUM_SAMPLES*OVERSAMPLE))); % Oversampled response function
    imp_resp = imp_resp/sum(output_meta.Grid.Row.WgtFunct); % Normalize to unit peak
    ind = find(imp_resp<1/sqrt(2),1,'first')+[-1 -0]; % Samples surrounding half-power point
    ind = interp1(imp_resp(ind), ind, 1/sqrt(2)); % Linear interpolation to solve for half-power point
    row_broadening_factor = 2*(ind - 1)/OVERSAMPLE;
    output_meta.Grid.Row.ImpRespWid = row_broadening_factor/output_meta.Grid.Row.ImpRespBW;
end
output_meta.Grid.Col.WgtType.WindowName = upper(char(xp.evaluate(...
    xpath_str({'product','imageGenerationParameters','sarProcessingInformation','azimuthWindow','windowName'}),...
    xml_domnode)));

%% Radar Collection
% Ultrafine and spotlight modes have "lower" and "upper" parts to the
% pulse.
% output_meta.RadarCollection.RefFreqIndex=uint32(0); % Absence of this field means all frequencies are true values
num_pulse_parts = str2double(xp.evaluate(...
    ['count(' xpath_str({'product','sourceAttributes','radarParameters','pulseLength'}) ')'],...
    xml_domnode));
for i=1:num_pulse_parts
    output_meta.RadarCollection.Waveform.WFParameters(i).TxRFBandwidth = str2double(xp.evaluate(...
        [xpath_str({'product','sourceAttributes','radarParameters','pulseBandwidth'}) '[' num2str(i) ']'],...
        xml_domnode)); % Bandwidth
    output_meta.RadarCollection.Waveform.WFParameters(i).TxPulseLength = str2double(xp.evaluate(...
        [xpath_str({'product','sourceAttributes','radarParameters','pulseLength'}) '[' num2str(i) ']'],...
        xml_domnode));
    output_meta.RadarCollection.Waveform.WFParameters(i).RcvDemodType='CHIRP';
    sample_rate = str2double(xp.evaluate(...
        [xpath_str({'product','sourceAttributes','radarParameters','adcSamplingRate'}) '[' num2str(i) ']'],...
        xml_domnode));
    output_meta.RadarCollection.Waveform.WFParameters(i).RcvWindowLength = str2double(xp.evaluate(...
        xpath_str({'product','sourceAttributes','radarParameters','samplesPerEchoLine'}),...
        xml_domnode))/sample_rate;
    output_meta.RadarCollection.Waveform.WFParameters(i).ADCSampleRate = sample_rate;
    output_meta.RadarCollection.Waveform.WFParameters(i).RcvFMRate = 0; % True for RcvDemodType='CHIRP'
end
bw = sum([output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth]);
output_meta.RadarCollection.TxFrequency.Min = fc-(bw/2); % fc calculated in Grid section
output_meta.RadarCollection.TxFrequency.Max = fc+(bw/2);
% Assumes pulse parts are exactly adjacent in bandwidth
output_meta.RadarCollection.Waveform.WFParameters(1).TxFreqStart = ...
    output_meta.RadarCollection.TxFrequency.Min;
for i=2:num_pulse_parts
    output_meta.RadarCollection.Waveform.WFParameters(i).TxFreqStart = ...
        output_meta.RadarCollection.Waveform.WFParameters(i-1).TxFreqStart + ...
        output_meta.RadarCollection.Waveform.WFParameters(i-1).TxRFBandwidth;
end
% Polarization
pols = textscan(char(xp.evaluate(...
    xpath_str({'product','sourceAttributes','radarParameters','polarizations'}),...
    xml_domnode)),'%s');
pols = pols{1};
tx_pols = unique(cellfun(@(x) x(1), pols));
for i=1:numel(pols)
    output_meta.RadarCollection.RcvChannels.ChanParameters(i).TxRcvPolarization = ...
        [pols{i}(1) ':' pols{i}(2)];
end
if isscalar(tx_pols) % Only one transmit polarization
    output_meta.RadarCollection.TxPolarization = tx_pols;
else % Multiple transmit polarizations
    output_meta.RadarCollection.TxPolarization = 'SEQUENCE';
    for i = 1:numel(tx_pols)
        output_meta.RadarCollection.TxSequence.TxStep(i).TxPolarization = tx_pols(i);
    end
end
% Another way to get polarimetric channels:
% num_pol_bands = str2double(xp.evaluate(...
%     ['count(' xpath_str({'product','imageAttributes','fullResolutionImageData'}) ')'],...
%     xml_domnode));
% for i=1:num_pol_bands
%     pol = char(xp.evaluate(...
%         [xpath_str({'product','imageAttributes','fullResolutionImageData'}) '[' num2str(i) ']/@pole'],...
%         xml_domnode));
% end

%% Timeline
output_meta.Timeline.CollectStart = rawDataStartTime + (rawDataStartTimeFrac/SECONDS_IN_A_DAY);
if strcmp(gen, 'RS2')
    prf_xp_str = {'pulseRepetitionFrequency'};
elseif strcmp(gen, 'RCM')
    prf_xp_str = {'prfInformation', 'pulseRepetitionFrequency'};
end
prf = str2double(xp.evaluate(...
    xpath_str([{'product','sourceAttributes', 'radarParameters'} prf_xp_str]),...
    xml_domnode));
num_lines_entries = str2double(xp.evaluate(...
    ['count(' xpath_str({'product','imageGenerationParameters','sarProcessingInformation','numberOfLinesProcessed'}) ')'],...
    xml_domnode));
num_lines_processed = zeros(num_lines_entries,1);
for i = 1:num_lines_entries
    num_lines_processed(i) = str2double(xp.evaluate(...
        [xpath_str({'product','imageGenerationParameters','sarProcessingInformation','numberOfLinesProcessed'}) '[@pole="' pols{i} '"]'],...
        xml_domnode));
end
num_lines_processed = num_lines_processed(1) * numel(tx_pols);
if num_lines_entries ~= numel(pols) || ~all(num_lines_processed == num_lines_processed(1))
    % This should never happen, but we'll throw an error if it does.
    warning('META2SICD_RS_XML:UnableToComputeCollectDuration', 'Unhandled data condition.');
end
prf = prf * num_pulse_parts;
if num_pulse_parts==2 && ...
        strcmp(output_meta.CollectionInfo.RadarMode.ModeType,'STRIPMAP')
    % Why????
    % This seems to be necessary to make CollectDuration match CA ranges
    % and to make 1/prf roughly equal to ss_zd_s (which is generally true
    % for STRIPMAP).  But we already doubled the prf above to account for
    % the pulse parts (to make real vs effective prf), so why do we have to
    % do it again? And why don't we have to do it for SPOTLIGHT?
    prf = 2*prf;
end
output_meta.Timeline.CollectDuration = num_lines_processed/prf;
output_meta.Timeline.IPP.Set.TStart = 0;
output_meta.Timeline.IPP.Set.TEnd = 0; % Apply real value later.  Just a placeholder.
output_meta.Timeline.IPP.Set.IPPStart = uint32(0);
output_meta.Timeline.IPP.Set.IPPEnd = uint32(num_lines_processed);
output_meta.Timeline.IPP.Set.IPPPoly = [0; prf];
output_meta.Timeline.IPP.Set.TEnd = output_meta.Timeline.CollectDuration;

%% Image Formation
output_meta.ImageFormation.RcvChanProc = ...
    struct('NumChanProc', uint32(1), ... % Assumes not a MODEX collect
    'PRFScaleFactor', 1/max(num_pulse_parts, numel(tx_pols))); % Either polarimetric or multi-step, but not both.
output_meta.ImageFormation.ImageFormAlgo = 'RMA';
output_meta.ImageFormation.TStartProc = 0;
output_meta.ImageFormation.TEndProc = output_meta.Timeline.CollectDuration;
output_meta.ImageFormation.TxFrequencyProc.MinProc = ...
    output_meta.RadarCollection.TxFrequency.Min;
output_meta.ImageFormation.TxFrequencyProc.MaxProc = ...
    output_meta.RadarCollection.TxFrequency.Max;
output_meta.ImageFormation.STBeamComp = 'NO';
output_meta.ImageFormation.ImageBeamComp = 'NO';
output_meta.ImageFormation.AzAutofocus = 'NO';
output_meta.ImageFormation.RgAutofocus = 'NO';

%% RMA.INCA
output_meta.RMA.RMAlgoType = 'OMEGA_K';
output_meta.RMA.ImageType = 'INCA';
output_meta.SCPCOA.SideOfTrack = char(xp.evaluate(...
    xpath_str({'product','sourceAttributes','radarParameters','antennaPointing'}),...
    xml_domnode));  % Should always be right looking for RCM
output_meta.SCPCOA.SideOfTrack = upper(output_meta.SCPCOA.SideOfTrack(1));
if output_meta.SCPCOA.SideOfTrack=='L'
    ss_zd_s = -ss_zd_s;
    % In addition to left/right, RS2 data can independently be in
    % increasing/decreasing line order.
    if (round((zd_first-zd_last)*SECONDS_IN_A_DAY) + ...
        (zd_first_frac-zd_last_frac)) < 0 % zd_last occurred after zd_first
        zd_first = zd_last;
        zd_first_frac = zd_last_frac;
    end
    look = 1;
else
    if (round((zd_first-zd_last)*SECONDS_IN_A_DAY) + ...
        (zd_first_frac-zd_last_frac)) > 0 % zd_last occurred before zd_first
        zd_first = zd_last;
        zd_first_frac = zd_last_frac;
    end
    look = -1;
end
% Zero doppler time of SCP relative to collect start
zd_t_scp = round((zd_first-rawDataStartTime)*SECONDS_IN_A_DAY) + ... % Convert days to seconds
    (zd_first_frac-rawDataStartTimeFrac) + ... % Handle fractional seconds
    (double(output_meta.ImageData.SCPPixel.Col) * ss_zd_s);
if strcmp(gen, 'RS2')
    near_range = str2double(xp.evaluate(...
        xpath_str({'product','imageGenerationParameters','sarProcessingInformation','slantRangeNearEdge'}),...
        xml_domnode)); % in meters
elseif strcmp(gen, 'RCM')
    near_range = str2double(xp.evaluate(...
        xpath_str({'product','sceneAttributes','imageAttributes','slantRangeNearEdge'}),...
        xml_domnode)); % in meters
end
output_meta.RMA.INCA.R_CA_SCP = near_range + ...
    (double(output_meta.ImageData.SCPPixel.Row)*output_meta.Grid.Row.SS);
output_meta.RMA.INCA.FreqZero = fc;
% Doppler Rate (We do this first since some other things are dependent on
% it.)
% For the purposes of the DRateSFPoly computation, we ignore any
% changes in velocity over the azimuth dimension.
pos_coefs = [P_x(:) P_y(:) P_z(:)];
% Velocity is derivate of position.
vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
vel_x = polyval(vel_coefs(:,1), zd_t_scp);
vel_y = polyval(vel_coefs(:,2), zd_t_scp);
vel_z = polyval(vel_coefs(:,3), zd_t_scp);
vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
r_ca = [output_meta.RMA.INCA.R_CA_SCP; 1]; % Polynomial representing range as a function of range distance from SCP
if strcmp(gen, 'RS2')
    drc_xp_str = {'product','imageGenerationParameters','dopplerRateValues', 'dopplerRateValuesCoefficients'};
elseif strcmp(gen, 'RCM')
    drc_xp_str = {'product','dopplerRate','dopplerRateEstimate', 'dopplerRateCoefficients'};
end
drr_xp_str = [drc_xp_str{1:(end-1)} {'dopplerRateReferenceTime'}];
dop_rate_coefs = str2num(xp.evaluate(xpath_str(drc_xp_str),... % Multiple numbers.  We need str2num instead of str2double
    xml_domnode)); % Shifted (origin at dop_rate_ref_t, not SCP) and scaled (sec, not m) version of SICD DopCentroidPoly
dop_rate_ref_t = str2double(xp.evaluate(xpath_str(drr_xp_str),...
    xml_domnode)); % Reference time of Doppler rate polynomial
dop_rate_coefs_shifted = polyshift(dop_rate_coefs, ... % Shift so SCP is reference
    (output_meta.RMA.INCA.R_CA_SCP*2/SPEED_OF_LIGHT) - ... % SICD reference time (SCP)
    dop_rate_ref_t); % Reference time of native Doppler Centroid polynomial
dop_rate_coefs_scaled = dop_rate_coefs_shifted .*  ... % Scale from seconds to meters
    (2/SPEED_OF_LIGHT) .^ (0:(length(dop_rate_coefs)-1));
output_meta.RMA.INCA.DRateSFPoly = - conv(dop_rate_coefs_scaled.',r_ca) * ... % Multiplication of two polynomials is just a convolution of their coefficients
    SPEED_OF_LIGHT / (2 * fc * vm_ca_sq(1)); % Assumes a SGN of -1

%% Fields dependent on Doppler rate
% This computation of SS is actually better than the claimed SS
% (sampledLineSpacing) in many ways, because this makes all of the metadata
% internally consistent.  This must be the sample spacing exactly at SCP
% (which is the definition for SS in SICD), if the other metadata from
% which is it computed is correct and consistent. Since column SS can vary
% slightly over a RGZERO image, we don't know if the claimed sample spacing
% in the native metadata is at our chosen SCP, or another point, or an
% average across image or something else.
output_meta.Grid.Col.SS = sqrt(vm_ca_sq(1)) * ss_zd_s * ...
    output_meta.RMA.INCA.DRateSFPoly(1,1);
output_meta.Grid.Col.ImpRespBW = dop_bw*abs(ss_zd_s)/output_meta.Grid.Col.SS; % Convert to azimuth spatial bandwidth (cycles per meter)
output_meta.RMA.INCA.TimeCAPoly = [zd_t_scp; ss_zd_s/output_meta.Grid.Col.SS];
if strcmpi(output_meta.Grid.Col.WgtType.WindowName,'KAISER') % The usual RS2 weigting
    output_meta.Grid.Col.WgtType.Parameter.name = 'BETA';
    output_meta.Grid.Col.WgtType.Parameter.value = char(xp.evaluate(...
        xpath_str({'product','imageGenerationParameters','sarProcessingInformation','azimuthWindow','windowCoefficient'}),...
        xml_domnode));
    beta_col = str2double(output_meta.Grid.Col.WgtType.Parameter.value);
    % We don't use the Mathworks Kaiser function, so we won't be dependent on the Signal Processing Toolbox
    output_meta.Grid.Col.WgtFunct = kaiser_nosptb(NUM_SAMPLES,beta_col);
    imp_resp = abs(fft(output_meta.Grid.Col.WgtFunct, round(NUM_SAMPLES*OVERSAMPLE))); % Oversampled response function
    imp_resp = imp_resp/sum(output_meta.Grid.Col.WgtFunct); % Normalize to unit peak
    ind = find(imp_resp<1/sqrt(2),1,'first')+[-1 -0]; % Samples surrounding half-power point
    ind = interp1(imp_resp(ind), ind, 1/sqrt(2)); % Linear interpolation to solve for half-power point
    col_broadening_factor = 2*(ind - 1)/OVERSAMPLE;
    output_meta.Grid.Col.ImpRespWid = col_broadening_factor/output_meta.Grid.Col.ImpRespBW;
end

%% Doppler Centroid
if strcmp(gen, 'RS2')
    dc_xp_str = {'product','imageGenerationParameters','dopplerCentroid'};
elseif strcmp(gen, 'RCM')
    dc_xp_str = {'product','dopplerCentroid','dopplerCentroidEstimate'};
end
dop_cent_coefs = str2num(xp.evaluate(... % Multiple numbers.  We need str2num instead of str2double
    xpath_str([dc_xp_str,{'dopplerCentroidCoefficients'}]),...
    xml_domnode)); % Shifted (origin at dop_cent_ref_t, not SCP) and scaled (sec, not m) version of SICD DopCentroidPoly
dop_cent_ref_t = str2double(xp.evaluate(...
    xpath_str([dc_xp_str,{'dopplerCentroidReferenceTime'}]),...
    xml_domnode)); % Reference time of Doppler Centroid polynomial
dop_cent_coefs_shifted = polyshift(dop_cent_coefs, ... % Shift so SCP is reference
    (output_meta.RMA.INCA.R_CA_SCP*2/SPEED_OF_LIGHT) - ... % SICD reference time (SCP)
    dop_cent_ref_t); % Reference time of native Doppler Centroid polynomial
dop_cent_coefs_scaled = dop_cent_coefs_shifted .*  ... % Scale from seconds to meters
    (2/SPEED_OF_LIGHT) .^ (0:(length(dop_cent_coefs)-1));
output_meta.RMA.INCA.DopCentroidPoly=dop_cent_coefs_scaled.';
% Adjust Doppler Centroid for spotlight
if strcmp(output_meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    [dop_est, dop_est_frac] = datenum_w_frac(char(xp.evaluate(...  % Doppler estimate time
        xpath_str([dc_xp_str 'timeOfDopplerCentroidEstimate']), xml_domnode)));
    dop_est_t = abs(round((dop_est - rawDataStartTime)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (dop_est_frac - rawDataStartTimeFrac)); % Handle fractional seconds
    dop_est_col = (dop_est_t - zd_t_scp)/ss_zd_s; % This is the column where the doppler centroid was computed.
    % Column-dependent variation in DopCentroidPoly due to spotlight
    output_meta.RMA.INCA.DopCentroidPoly(1,2) = ...
        -look * fc * (2 / SPEED_OF_LIGHT) * sqrt(vm_ca_sq(1)) / ...
        output_meta.RMA.INCA.R_CA_SCP;
    % dopplerCentroid in native metadata was defined at specific column,
    % which might not be our SCP column.  Adjust so that SCP column is
    % correct.
    output_meta.RMA.INCA.DopCentroidPoly(1,1) = ...
        output_meta.RMA.INCA.DopCentroidPoly(1,1) - ...
        (output_meta.RMA.INCA.DopCentroidPoly(1,2) * ...
        dop_est_col * output_meta.Grid.Col.SS);
end
output_meta.Grid.Col.DeltaKCOAPoly = ...
    output_meta.RMA.INCA.DopCentroidPoly * ss_zd_s / output_meta.Grid.Col.SS;
% Compute Col.DeltaK1/K2 from DeltaKCOAPoly
% This is not always straightforward to do generically for any possible
% DeltaKCOAPoly 2D polynomial, since you would have to compute all 2D roots
% and edge cases.  However, for the RS case, this can be solved exactly,
% since its usually a 1D polynomial.  Even the spotlight case only brings
% in a linear variation in the second dimensions, so its still easily
% solved.
% Min/max in row/range must exist at edges or internal local min/max
minmax = roots(polyder(output_meta.Grid.Col.DeltaKCOAPoly(end:-1:1,1)));
rg_bounds_m = (double([0 (output_meta.ImageData.NumRows-1)]) - ...
     double(output_meta.ImageData.SCPPixel.Row)) * output_meta.Grid.Row.SS;
possible_bounds_rg = [rg_bounds_m minmax(minmax>min(rg_bounds_m) & minmax<max(rg_bounds_m))];
% Constant or (linearly increasing\decreasing for spotlight) in column, so
% edges must contain max/min.
possible_bounds_az = (double([0 (output_meta.ImageData.NumCols-1)]) - ...
     double(output_meta.ImageData.SCPPixel.Col)) * output_meta.Grid.Col.SS;
possible_bounds_deltak = sicd_polyval2d(output_meta.Grid.Col.DeltaKCOAPoly, ...
    possible_bounds_az, possible_bounds_rg);
output_meta.Grid.Col.DeltaK1 = min(possible_bounds_deltak(:)) - ...
    (output_meta.Grid.Col.ImpRespBW/2);
output_meta.Grid.Col.DeltaK2 = max(possible_bounds_deltak(:)) + ...
    (output_meta.Grid.Col.ImpRespBW/2);
% Wrapped spectrum
if (output_meta.Grid.Col.DeltaK1 < -(1/output_meta.Grid.Col.SS)/2) || ...
        (output_meta.Grid.Col.DeltaK2 > (1/output_meta.Grid.Col.SS)/2)
    output_meta.Grid.Col.DeltaK1 = -(1/output_meta.Grid.Col.SS)/2;
    output_meta.Grid.Col.DeltaK2 = -output_meta.Grid.Col.DeltaK1;
end
%% TimeCOAPoly
% TimeCOAPoly=TimeCA+(DopCentroid/dop_rate)
% Since we can't evaluate this equation analytically, we will evaluate
% samples of it across our image and fit a 2D polynomial to it.
POLY_ORDER = 2; % Order of polynomial which we want to compute
grid_samples = POLY_ORDER + 1; % in each dimension
coords_az_m = linspace(-double(output_meta.ImageData.SCPPixel.Col) * output_meta.Grid.Col.SS,...
    double(output_meta.ImageData.NumCols-output_meta.ImageData.SCPPixel.Col) * ...
    output_meta.Grid.Col.SS,grid_samples);
coords_rg_m = linspace(-double(output_meta.ImageData.SCPPixel.Row) * output_meta.Grid.Row.SS,...
    double(output_meta.ImageData.NumRows-output_meta.ImageData.SCPPixel.Row) * ...
    output_meta.Grid.Row.SS, grid_samples);
timeca_sampled = sicd_polyval2d(output_meta.RMA.INCA.TimeCAPoly(:).',coords_az_m,coords_rg_m);
dopcentroid_sampled = sicd_polyval2d(output_meta.RMA.INCA.DopCentroidPoly,coords_az_m,coords_rg_m);
doprate_sampled = sicd_polyval2d(dop_rate_coefs_scaled.',coords_az_m,coords_rg_m);
timecoapoly_sampled = timeca_sampled+(dopcentroid_sampled./doprate_sampled);
% Least squares fit for 2D polynomial
% A*x = b
[coords_az_m, coords_rg_m] = ndgrid(coords_az_m, coords_rg_m);
a = zeros(grid_samples^2, (POLY_ORDER+1)^2);
for i = 0:POLY_ORDER
    for j = 0:POLY_ORDER
        a(:,i*(POLY_ORDER+1)+j+1) = (coords_rg_m(:).^j).*(coords_az_m(:).^i);
    end
end
b_coa = zeros((POLY_ORDER+1)^2,1);
for i=1:((POLY_ORDER+1)^2)
   b_coa(i)=sum(timecoapoly_sampled(:).*a(:,i)); % center of aperture
end
A=zeros((POLY_ORDER+1)^2);
for i=1:((POLY_ORDER+1)^2)
    for j=1:((POLY_ORDER+1)^2)
        A(i,j)=sum(a(:,i).*a(:,j));
    end
end
old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
x=A\b_coa; % MATLAB often flags this as badly scaled, but results still appear valid
warning(old_warning_state);
output_meta.Grid.TimeCOAPoly=reshape(x, POLY_ORDER+1, POLY_ORDER+1);
if strcmp(output_meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    output_meta.Grid.TimeCOAPoly = output_meta.Grid.TimeCOAPoly(1);
    % This field required to compute TimeCOAPoly, but not allowed for
    % spotlight in SICD.
    output_meta.RMA.INCA = rmfield(output_meta.RMA.INCA, 'DopCentroidPoly');
else % This field also not allowed for spotlight in SICD.
    output_meta.RMA.INCA.DopCentroidCOA=true;
end

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

%% Radiometric
if exist('beta_domnode','var')
    % Offset always zero for SLC, so we only need gains
    betas = str2num(xp.evaluate(xpath_str({'lut','gains'}),beta_domnode)); %#ok<ST2NM>
    % Of the provided LUTs, we really only work with beta here, since it is
    % the radiometric term most commonly kept constant, and the others
    % (sigma/gamma) can always be derived from it.
    % if any(strcmp({'Constant-beta', 'Point Target', 'Point Target-1', ...
    %         'Calibration-1', 'Calibration-2', 'Ship-1', 'Ship-2', ...
    %         'Ship-3', 'Unity'}, ... % Known modes with constant beta
    %         xp.evaluate(xpath_str({'product','imageGenerationParameters', ...
    %         'sarProcessingInformation','lutApplied'}), xml_domnode)))
    if all(betas==betas(1)) % In case we missed some modes, this condition may be more reliable
        output_meta.Radiometric.BetaZeroSFPoly = 1/betas(1)^2;
    else  % Otherwise fit a 1D polynomial in range
        % RS2 has value for every row
        coords_rg_m = (double(0:(output_meta.ImageData.NumRows-1)) - ...
            double(output_meta.ImageData.SCPPixel.Row)) * output_meta.Grid.Row.SS;
        if strcmp(gen, 'RCM') % RCM subsamples the rows
            rng_indices = ((str2double(xp.evaluate(xpath_str({'lut','pixelFirstAnglesValue'}), beta_domnode)):...
                ... % Should be this, but simulated datasets are inconsistent with spec document
                ... % crm_indices = (xp.evaluate(xpath_str({'lut','pixelFirstLutValue'}), beta_domnode):...
                (str2double(xp.evaluate(xpath_str({'lut','numberOfValues'}), beta_domnode))-1)) * ... % First row is zero
                str2double(xp.evaluate(xpath_str({'lut','stepSize'}), beta_domnode))) + 1;
            coords_rg_m = coords_rg_m(rng_indices);
        end
        % For the datasets we have seen, this function is very close to
        % linear.  For the "Mixed" LUT, there is a tiny linear piecewise
        % deviation from a single overall linear.
        betapoly = polyfit(coords_rg_m, betas, 2);
        output_meta.Radiometric.BetaZeroSFPoly = betapoly(end:-1:1).';
    end
    % RCS, Sigma, and Gamma will be computed below in derived_sicd_fields
    % derived_sicd_fields.
    output_meta.Radiometric.NoiseLevel.NoiseLevelType = 'ABSOLUTE';
    if strcmp(gen, 'RS2')
        beta0_str = [xpath_str({'product', 'sourceAttributes', ...
            'radarParameters','referenceNoiseLevel'}) ...
            '[@incidenceAngleCorrection="Beta Nought"]'];
        noise_domnode = xml_domnode;  % RS2 noise is in main product.xml
    elseif strcmp(gen, 'RCM') % RCM noise is in separate file
        beta0_str = [xpath_str({'noiseLevels', 'referenceNoiseLevel'}) ...
            '/*[text()=''Beta Nought'']/..'];
    end
    if exist('noise_domnode','var')
        pfv = str2double(xp.evaluate([beta0_str ...
            '/*[local-name()=''pixelFirstNoiseValue'']'], noise_domnode)); % Index of first row defined to be zero
        step = str2double(xp.evaluate([beta0_str ...
            '/*[local-name()=''stepSize'']'], noise_domnode));
        beta0s = str2num(xp.evaluate([beta0_str ...
            '/*[local-name()=''noiseLevelValues'']'], noise_domnode)); %#ok<ST2NM>
        range_coords = output_meta.Grid.Row.SS * ...
            ((((1:numel(beta0s))-1) * step) + pfv - double(output_meta.ImageData.SCPPixel.Row));
        noisepoly = polyfit(range_coords, beta0s - (10*log10(polyval(...
            output_meta.Radiometric.BetaZeroSFPoly(end:-1:1), range_coords))), 2);
        output_meta.Radiometric.NoiseLevel.NoisePoly = noisepoly(end:-1:1).';
    end
end

%% SCPCOA
% All of these fields are derivable for more fundamental fields.
output_meta = derived_sicd_fields(output_meta);

%% Process fields specific to each polarimetric band
band_independent_meta = output_meta; % Values that are consistent across all bands
grouped_meta = cell(numel(pols),1);
for i=1:numel(pols)
    output_meta = band_independent_meta;
    
    output_meta.ImageFormation.RcvChanProc.ChanIndex = i;
    output_meta.ImageFormation.TxRcvPolarizationProc = ...
        output_meta.RadarCollection.RcvChannels.ChanParameters(i).TxRcvPolarization;
    
    grouped_meta{i} = output_meta;
end
output_meta = grouped_meta; % Cell array with metadata struct for each band

end

% Creates an XPath query that is a namespace insensitive search for an XML
% heirachy.  Input is a cell array of strings specifying the XML fields,
% starting with the root node and working its way down the XML tree to the
% leaf nodes.  This is required since RS2 XML files specify a default
% namespace, so XPath queries with no namespace speficied will not work.
% (That only searches for elements associated with no namespace at all.)
% Note: This seems to be an issue only in MATLAB 2010a and above, as
% previous versions used another XPATH library that was more tolerant.
function out_str = xpath_str(in_cell_array)
    out_str='';
    for j=1:length(in_cell_array)
        out_str=[out_str '/*[local-name()=''' in_cell_array{j} ''']'];
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
    if isnan(datenum_frac), datenum_frac = 0; end;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////