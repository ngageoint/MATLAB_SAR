function [ output_meta ] = meta2sicd_s1product( domnode, eof_domnode )
%META2SICD_S1PRODUCT Converts Sentinel-1 annotation "product" XML file into SICD format
%
% Written by: Wade Schwartzkopf, NGA Research
%
% We do not populate the SICD Antenna field.  An antenna pattern is given
% in the metadata, but its not clear what you would use this for, so we
% haven't bothered to convert to the SICD structure yet.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup
SECONDS_IN_A_DAY = 24*60*60;
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

%% CollectionInfo
common_meta.CollectionInfo.CollectorName = char(xp.evaluate(...
    'product/adsHeader/missionId',domnode));
common_meta.CollectionInfo.CollectType='MONOSTATIC';
common_meta.CollectionInfo.RadarMode.ModeID = char(xp.evaluate(...
    'product/adsHeader/mode',domnode));
if common_meta.CollectionInfo.RadarMode.ModeID(1)=='S'
    common_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
else
    % Actually TOPSAR.  Not what we normally think of for DYNAMIC STRIPMAP,
    % but it is definitely not SPOTLIGHT (actually counter to the spotlight
    % beam motion), and it isn't STRIPMAP with a constant angle between the
    % beam and direction of travel either, so we use DYNAMIC STRIPMAP as a
    % catchall.
    common_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
end
common_meta.CollectionInfo.Classification='UNCLASSIFIED';

%% ImageData
% For SLC, the following test should always hold true:
if strcmpi(xp.evaluate('product/imageAnnotation/imageInformation/pixelValue',domnode),'Complex') && ...
    strcmpi(xp.evaluate('product/imageAnnotation/imageInformation/outputPixels',domnode),'16 bit Signed Integer')
    common_meta.ImageData.PixelType = 'RE16I_IM16I';
else
    warning('META2SICD_S1PRODUCT:UNRECOGNIZED_SLC','SLC data should be 16-bit complex.');
end
num_bursts = str2double(xp.evaluate(...
    'count(product/swathTiming/burstList/burst)',domnode));
% These two definitions of NumRows should always be the same for
% non-STRIPMAP data (For STRIPMAP, samplesPerBurst is set to zero.)  Number
% of rows in burst should be the same as the full image.  Both of these
% numbers also should match the ImageWidth field of the measurement TIFF.
if num_bursts>0
    common_meta.ImageData.NumRows = uint32(str2double(char(xp.evaluate(...
        'product/swathTiming/samplesPerBurst',domnode))));
else % STRIPMAP
    common_meta.ImageData.NumRows = uint32(str2double(xp.evaluate(...
        'product/imageAnnotation/imageInformation/numberOfSamples',domnode)));
end
% These NumCols definition will be different.  Since each burst is its own
% coherent data period, and thus SICD, we set the SICD metadata to describe
% each individual burst.
if num_bursts>0
    % Ths in the number of coumns in a single burst.
    common_meta.ImageData.NumCols = uint32(str2double(char(xp.evaluate(...
        'product/swathTiming/linesPerBurst',domnode))));
else % STRIPMAP
    % This in the number of columns in the full TIFF measurement file, even
    % if it contains multiple bursts.
    common_meta.ImageData.NumCols = uint32(str2double(xp.evaluate(...
        'product/imageAnnotation/imageInformation/numberOfLines',domnode)));
end

common_meta.ImageData.FirstRow = uint32(0);
common_meta.ImageData.FirstCol = uint32(0);
common_meta.ImageData.FullImage.NumRows = common_meta.ImageData.NumRows;
common_meta.ImageData.FullImage.NumCols = common_meta.ImageData.NumCols;
% SCP pixel within entire TIFF
center_cols = round((-0.5 + (1:max(num_bursts,1))) * double(common_meta.ImageData.NumCols)) - 1;
center_rows = repmat(round(double(common_meta.ImageData.NumRows)/2 - 1),size(center_cols));
% SCP pixel within single burst image is the same for all burst, since east
% burst is the same size
common_meta.ImageData.SCPPixel.Col = center_cols(1);
common_meta.ImageData.SCPPixel.Row = center_rows(1);

%% GeoData
common_meta.GeoData.EarthModel = 'WGS_84';
% Initially, we just seed this with a rough value.  Later we will put
% in something more precise.
num_grid_points=str2double(xp.evaluate(...
    'count(product/geolocationGrid/geolocationGridPointList/geolocationGridPoint)',domnode));
[scp_col, scp_row, lat, lon, hgt] = deal(zeros(num_grid_points,1));
for j = 1:num_grid_points
    scp_col(j) = str2double(xp.evaluate(...
        ['product/geolocationGrid/geolocationGridPointList/geolocationGridPoint[' num2str(j) ']/line'],...
        domnode));
    scp_row(j) = str2double(xp.evaluate(...
        ['product/geolocationGrid/geolocationGridPointList/geolocationGridPoint[' num2str(j) ']/pixel'],...
        domnode));
    lat(j) = str2double(xp.evaluate(...
        ['product/geolocationGrid/geolocationGridPointList/geolocationGridPoint[' num2str(j) ']/latitude'],...
        domnode));
    lon(j) = str2double(xp.evaluate(...
        ['product/geolocationGrid/geolocationGridPointList/geolocationGridPoint[' num2str(j) ']/longitude'],...
        domnode));
    hgt(j) = str2double(xp.evaluate(...
        ['product/geolocationGrid/geolocationGridPointList/geolocationGridPoint[' num2str(j) ']/height'],...
        domnode));
end
scp_lats = griddata(scp_col, scp_row, lat, center_cols, center_rows);
scp_lons = griddata(scp_col, scp_row, lon, center_cols, center_rows);
scp_hgts = griddata(scp_col, scp_row, hgt, center_cols, center_rows);

%% Grid
common_meta.Grid.ImagePlane = upper(strtok(char(xp.evaluate(...
    'product/generalAnnotation/productInformation/projection',domnode))));
common_meta.Grid.Type = 'RGZERO';
delta_tau_s = 1/str2double(char(xp.evaluate(...
    'product/generalAnnotation/productInformation/rangeSamplingRate',domnode)));
common_meta.Grid.Row.SS = (SPEED_OF_LIGHT/2) * delta_tau_s;
% common_meta.Grid.Row.SS = str2double(xp.evaluate(... % Another way to do it
%     'product/imageAnnotation/imageInformation/rangePixelSpacing',domnode));
common_meta.Grid.Row.Sgn = -1;
% Justification for Sgn:
% 1) "Sentinel-1 Level 1 Detailed Algorithm Definition" shows last step in
% image formation as IFFT, which would mean a forward FFT (-1 Sgn) would be
% required to transform back.
% 2) The forward FFT of a sliding window shows the Doppler centroid
% increasing as you move right in the image, which must be the case for the
% TOPSAR collection mode which starts in a rear squint and transitions to a
% forward squint (and are always right looking).
fc = str2double(char(xp.evaluate(...
    'product/generalAnnotation/productInformation/radarFrequency',domnode)));
common_meta.Grid.Row.KCtr = 2*fc/SPEED_OF_LIGHT;
common_meta.Grid.Row.DeltaKCOAPoly=0;
spp_str = 'product/imageAnnotation/processingInformation/swathProcParamsList/swathProcParams/';
common_meta.Grid.Row.ImpRespBW=2*str2double(xp.evaluate([spp_str 'rangeProcessing/processingBandwidth'],domnode))/SPEED_OF_LIGHT;
common_meta.Grid.Row.WgtType.WindowName=upper(char(xp.evaluate(...
    [spp_str 'rangeProcessing/windowType'],domnode)));
if strcmpi(common_meta.Grid.Row.WgtType.WindowName,'NONE')
    common_meta.Grid.Row.WgtType.WindowName = 'UNIFORM';
elseif strcmpi(common_meta.Grid.Row.WgtType.WindowName,'HAMMING') % The usual Sentinel weighting
    common_meta.Grid.Row.WgtType.Parameter.name = 'COEFFICIENT';
    common_meta.Grid.Row.WgtType.Parameter.value = char(xp.evaluate(...
        [spp_str 'rangeProcessing/windowCoefficient'],domnode));
    a = str2double(common_meta.Grid.Row.WgtType.Parameter.value); % Generalized Hamming window parameter
    common_meta.Grid.Row.WgtFunct = raised_cos_fun(512,a);
    % Computation of broadening factor for uniform window is:
    % 2 * fzero(@(x) (sin(pi*x)/(pi*x)) - (1/sqrt(2)), .1)
    row_broadening_factor = 2*fzero(@(x) ...
        a*(sin(pi*x)/(pi*x)) + ((1-a)*(sin(pi*(x-1))/(pi*(x-1)))/2) + ...
        ((1-a)*(sin(pi*(x+1))/(pi*(x+1)))/2) - a/sqrt(2), .1);
    common_meta.Grid.Row.ImpRespWid = row_broadening_factor/common_meta.Grid.Row.ImpRespBW;
end
common_meta.Grid.Col.SS = str2double(xp.evaluate(...
    'product/imageAnnotation/imageInformation/azimuthPixelSpacing',domnode));
common_meta.Grid.Col.Sgn = -1; % Must be the same as Row.Sgn
common_meta.Grid.Col.KCtr = 0;
dop_bw=str2double(xp.evaluate([spp_str 'azimuthProcessing/processingBandwidth'],domnode)); % Doppler bandwidth
ss_zd_s=str2double(xp.evaluate(... % Image column spacing in zero doppler time (seconds)
    'product/imageAnnotation/imageInformation/azimuthTimeInterval',...
    domnode)); % Sentinel-1 is always right-looking, so should always be positive
common_meta.Grid.Col.ImpRespBW=dop_bw*ss_zd_s/common_meta.Grid.Col.SS; % Convert to azimuth spatial bandwidth (cycles per meter)
common_meta.Grid.Col.WgtType.WindowName=upper(char(xp.evaluate(...
    [spp_str 'azimuthProcessing/windowType'],domnode)));
if strcmpi(common_meta.Grid.Col.WgtType.WindowName,'NONE')
    common_meta.Grid.Col.WgtType.WindowName = 'UNIFORM';
elseif strcmpi(common_meta.Grid.Col.WgtType.WindowName,'HAMMING') % The usual Sentinel weighting
    common_meta.Grid.Col.WgtType.Parameter.name = 'COEFFICIENT';
    common_meta.Grid.Col.WgtType.Parameter.value = char(xp.evaluate(...
        [spp_str 'azimuthProcessing/windowCoefficient'],domnode));
    a = str2double(common_meta.Grid.Col.WgtType.Parameter.value); % Generalized Hamming window parameter
    common_meta.Grid.Col.WgtFunct = raised_cos_fun(512,a);
    col_broadening_factor = 2*fzero(@(x) ...
        a*(sin(pi*x)/(pi*x)) + ((1-a)*(sin(pi*(x-1))/(pi*(x-1)))/2) + ...
        ((1-a)*(sin(pi*(x+1))/(pi*(x+1)))/2) - a/sqrt(2), .1);
    common_meta.Grid.Col.ImpRespWid = col_broadening_factor/common_meta.Grid.Col.ImpRespBW;
end
% We will compute Grid.Col.DeltaKCOAPoly separately per-burst later.
% The following Grid fields will be computed later in derived_sicd_fields
% Grid.Row/Col.WgtFunct
% Grid.Row/Col.ImpRespWid
% Grid.Row/Col.DeltaK1
% Grid.Row/Col.DeltaK2

%% Timeline
prf=str2double(xp.evaluate(...
    'product/generalAnnotation/downlinkInformationList/downlinkInformation/prf',domnode));
% Because of calibration pulses, it is unlikely this PRF was maintained
% through this entire period, but we don't currently include that detail.
common_meta.Timeline.IPP.Set.IPPPoly=[0; prf];
[azimuth_time_first_line, azimuth_time_frac] = datenum_w_frac(char(xp.evaluate(...
    'product/imageAnnotation/imageInformation/productFirstLineUtcTime',...
    domnode))); % Always the left-most SICD column (of first bursts or entire STRIPMAP dataset), since Sentinel-1 is always right-looking
eta_mid = ss_zd_s * double(common_meta.ImageData.SCPPixel.Col); % Offset in zero Doppler time from first column to SCP column

%% Position
% Compute state vectors
num_state_vectors=str2double(xp.evaluate(...
    'count(product/generalAnnotation/orbitList/orbit)',domnode));
state_vector_T  = zeros(num_state_vectors,1);
state_vector_T_frac  = zeros(num_state_vectors,1);
state_vector_pos  = zeros(num_state_vectors,3);
state_vector_vel = zeros(num_state_vectors,3);
for i=1:num_state_vectors
    timeStamp = char(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/time'],...
        domnode));
    [state_vector_T(i), state_vector_T_frac(i)] = datenum_w_frac(timeStamp);

    state_vector_pos(i,1) = str2double(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/position/x'],...
        domnode));
    state_vector_pos(i,2) = str2double(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/position/y'],...
        domnode));
    state_vector_pos(i,3) = str2double(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/position/z'],...
        domnode));
    state_vector_vel(i,1) = str2double(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/velocity/x'],...
        domnode));
    state_vector_vel(i,2) = str2double(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/velocity/y'],...
        domnode));
    state_vector_vel(i,3) = str2double(xp.evaluate(...
        ['product/generalAnnotation/orbitList/orbit[' num2str(i) ']/velocity/z'],...
        domnode));
end
% If external orbit file was passed in, use it instead of orbit state files
% in SLC annotation file
common_meta.CollectionInfo.Parameter{4}.name = 'ORBIT_SOURCE';
if exist('eof_domnode','var')
    % Use same time range as was used in SLC annotation file plus some
    % buffer to assure enough state vectors are selected, since we are
    % fitting to 5th order polynomial.  EOF files have state vectors every
    % 10 seconds.
    buffer = max(5, (70-(max(state_vector_T)-min(state_vector_T))*SECONDS_IN_A_DAY)/2);
    timerange = [min(state_vector_T - (buffer/SECONDS_IN_A_DAY)) ...
        max(state_vector_T + (buffer/SECONDS_IN_A_DAY))];
    [state_vector_T, state_vector_T_frac, state_vector_pos, state_vector_vel] = ...
        get_osv_from_eof(eof_domnode, timerange);
    
    common_meta.CollectionInfo.Parameter{4}.value = char(xp.evaluate(... % Orbit file type
        'Earth_Explorer_File/Earth_Explorer_Header/Fixed_Header/File_Type',eof_domnode));
else
    common_meta.CollectionInfo.Parameter{4}.value = 'SLC_INTERNAL';
end

%% RadarCollection
pol = char(xp.evaluate('product/adsHeader/polarisation',domnode));
common_meta.RadarCollection.TxPolarization = pol(1);
common_meta.RadarCollection.TxFrequency.Min = fc + str2double(xp.evaluate(...
    'product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/txPulseStartFrequency',domnode));
common_meta.RadarCollection.Waveform.WFParameters.TxFreqStart = ...
    common_meta.RadarCollection.TxFrequency.Min;
common_meta.RadarCollection.Waveform.WFParameters.TxPulseLength=str2double(xp.evaluate(...
    'product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/txPulseLength',domnode));
common_meta.RadarCollection.Waveform.WFParameters.TxFMRate=str2double(xp.evaluate(...
    'product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/txPulseRampRate',domnode));
bw = common_meta.RadarCollection.Waveform.WFParameters.TxPulseLength*...
    common_meta.RadarCollection.Waveform.WFParameters.TxFMRate;
common_meta.RadarCollection.TxFrequency.Max = common_meta.RadarCollection.TxFrequency.Min + bw;
common_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth=bw;
common_meta.RadarCollection.Waveform.WFParameters.RcvDemodType='CHIRP';
common_meta.RadarCollection.Waveform.WFParameters.RcvFMRate=0; % True for RcvDemodType='CHIRP'
common_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate=str2double(char(xp.evaluate(...
    'product/generalAnnotation/productInformation/rangeSamplingRate',domnode))); % Raw not decimated
% After decimation would be:
% output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate=str2double(xp.evaluate(...
%         'product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/rangeDecimation/samplingFrequencyAfterDecimation',domnode));
num_swl = str2double(xp.evaluate(...
        'count(product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/swlList/swl)',domnode));
for i = 1:num_swl % We could have multiple receive window lengths across the collect
    common_meta.RadarCollection.Waveform.WFParameters(i) = common_meta.RadarCollection.Waveform.WFParameters(1);
    common_meta.RadarCollection.Waveform.WFParameters(i).RcvWindowLength=str2double(xp.evaluate(...
        ['product/generalAnnotation/downlinkInformationList/downlinkInformation/downlinkValues/swlList/swl[' num2str(i) ']/value'],domnode));
end

%% ImageFormation
common_meta.ImageFormation.RcvChanProc=struct('NumChanProc',1,'PRFScaleFactor',1);
common_meta.ImageFormation.TxRcvPolarizationProc=[pol(1) ':' pol(2)];
% RcvChanProc.ChanIndex must be populated external to this since it depends
% on how the polarization were ordered in manifest file.
% Assume image formation uses all data
common_meta.ImageFormation.TStartProc=0;
common_meta.ImageFormation.TxFrequencyProc.MinProc=...
    common_meta.RadarCollection.TxFrequency.Min;
common_meta.ImageFormation.TxFrequencyProc.MaxProc=...
    common_meta.RadarCollection.TxFrequency.Max;
common_meta.ImageFormation.ImageFormAlgo='RMA';
% From the Sentinel-1 Level 1 Detailed Algorithm Definition document
if common_meta.CollectionInfo.RadarMode.ModeID(1)=='S'
    common_meta.ImageFormation.STBeamComp='NO'; %STRIPMAP
else
    common_meta.ImageFormation.STBeamComp='SV'; %TOPSAR
end
common_meta.ImageFormation.ImageBeamComp='NO';
common_meta.ImageFormation.AzAutofocus='NO';
common_meta.ImageFormation.RgAutofocus='NO';

%% RMA
% "Sentinel-1 Level 1 Detailed Algorithm Definition" document seems to most
% closely match the RangeDoppler algorithm (with accurate secondary range
% compression or "option 2" as described in the Cumming and Wong book).
common_meta.RMA.RMAlgoType = 'RG_DOP';
common_meta.RMA.ImageType = 'INCA';
tau_0 = str2double(char(xp.evaluate(... % tau_0 is notation from ESA deramping paper
    'product/imageAnnotation/imageInformation/slantRangeTime',domnode)));
common_meta.RMA.INCA.R_CA_SCP = (SPEED_OF_LIGHT/2) * (tau_0 + ...
        (double(common_meta.ImageData.SCPPixel.Row) * delta_tau_s));
common_meta.RMA.INCA.FreqZero = fc;
% If we use the Doppler Centroid as defined directly in the manifest.safe
% metadata, then the center of frequency support Col.DeltaKCOAPoly does not
% correspond to RMA.INCA.DopCentroidPoly.  However, we will compute
% TimeCOAPoly later to match a newly computed Doppler Centroid based off of
% DeltaKCOAPoly, assuming that the the COA is at the peak signal (fdop_COA
% = fdop_DC).
common_meta.RMA.INCA.DopCentroidCOA = true;

%% Doppler centroid
% Get common (non-burst specific) parameters we will need for Doppler
% centroid and rate computations later
num_dc_estimates_az=str2double(xp.evaluate(...
    'count(product/dopplerCentroid/dcEstimateList/dcEstimate)',domnode));
for i = 1:num_dc_estimates_az
    [dc_az_time_s(i), dc_az_time_frac(i)] = datenum_w_frac(char(xp.evaluate(...
        ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/azimuthTime'],...
        domnode)));
    % How to interpret the various representations of Doppler Centroid in
    % the Sentinel 1 XML:
    % This implementation allows for estimates at different azimuth times
    % to have different numbers of fineDce values, but I don't know if that
    % ever actually happens.
    % num_dc_estimates_rg=str2double(xp.evaluate(...
    %     ['count(product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/fineDceList/fineDce)'],domnode));
    % [tSR, freq] = deal(zeros(num_dc_estimates_rg,1));
    % for j = 1:num_dc_estimates_rg
    %     tSR(j) = str2double(xp.evaluate(... % Two-way slant range (s)
    %         ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/fineDceList/fineDce[' num2str(j) ']/slantRangeTime'],...
    %         domnode));
    %     freq(j) = str2double(xp.evaluate(... % Doppler estimate
    %         ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/fineDceList/fineDce[' num2str(j) ']/frequency'],...
    %         domnode));
    % end
    dc_t0(i) = str2double(xp.evaluate(...
        ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/t0'],...
        domnode));
    % % data_dc_poly should be equal to:
    % % polyfit(tSR-dc_t0, freq, 2)
    data_dc_poly{i} = str2num(xp.evaluate(...
        ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/dataDcPolynomial'],...
        domnode));
    % % rms_dc_poly should be equal to:
    % % rms(freq - polyval(polyfit(tSR-dc_t0, freq, 2), tSR-dc_t0))
    % rms_dc_poly = str2num(xp.evaluate(...
    %     ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/dataDcRmsError'],...
    %     domnode));
    % % An alternative method for computing the Doppler polynomial.  Uses the
    % % Doppler centroid estimate from orbit geometry.  Much smoother.
    % geom_dc_poly = str2num(xp.evaluate(...
    %     ['product/dopplerCentroid/dcEstimateList/dcEstimate[' num2str(i) ']/geometryDcPolynomial'],...
    %     domnode));
    % freq=polyval(geom_dc_poly(end:-1:1),tSR-dc_t0); % Samples at the same points as fineDceList
    
    % Of all these options, we will use dataDcPolynomial, since that's what
    % they use in the ESA document on TOPS SLC deramping.  The data-derived
    % samples are not always particularly smooth, nor are they sampled very
    % densely, so its not clear how accurate this polynomial representation
    % is.
end
num_azfm_rates=str2double(xp.evaluate(...
    'count(product/generalAnnotation/azimuthFmRateList/azimuthFmRate)',...
    domnode));
for i = 1:num_azfm_rates
    [az_t_s(i), az_t_frac(i)] = datenum_w_frac(char(xp.evaluate(...
        ['product/generalAnnotation/azimuthFmRateList/azimuthFmRate[' ...
        num2str(i) ']/azimuthTime'],domnode)));
    az_t0(i) = str2double(xp.evaluate(...
        ['product/generalAnnotation/azimuthFmRateList/azimuthFmRate[' ...
        num2str(i) ']/t0'],domnode)).';
    % Two different ways we have seen in XML for storing the FM Rate polynomial
    if strcmpi('true',xp.evaluate(...
            ['boolean(product/generalAnnotation/azimuthFmRateList/azimuthFmRate[' ...
            num2str(i) ']/azimuthFmRatePolynomial)'],domnode))
        k_a_poly{i} = str2num(xp.evaluate(...
            ['product/generalAnnotation/azimuthFmRateList/azimuthFmRate[' ...
            num2str(i) ']/azimuthFmRatePolynomial'],domnode)).';
    else
        for j=1:3 % Assume hard coded to 3
            k_a_poly{i}(j)=str2double(xp.evaluate(... % Doppler FM rate with regard to raw time
                ['product/generalAnnotation/azimuthFmRateList/azimuthFmRate[' ...
                num2str(i) ']/c' num2str(j-1) ],...
                domnode));
        end
    end
end

% Azimuth steering rate (constant, not dependent on burst or range)
k_psi = str2double(char(xp.evaluate(...
    'product/generalAnnotation/productInformation/azimuthSteeringRate',domnode)));
k_psi = k_psi*pi/180; % Convert from degrees/sec into radians/sec


%% Compute per/burst metadata
output_meta = cell(num_bursts,1);
for i = 1:max(num_bursts,1)
    output_meta{i} = common_meta;
    
    %% CollectionInfo
    % These values are needed to generate CoreName later, but we will be
    % have to wait until after Timeline to fully generate that string,
    % since we include date in it.
    if strcmpi('true',xp.evaluate('boolean(product/imageAnnotation/imageInformation/sliceNumber)',domnode))
        slice = str2double(xp.evaluate(...
            'product/imageAnnotation/imageInformation/sliceNumber',domnode));
    else
        slice = 0;
    end
    swath = char(xp.evaluate('product/adsHeader/swath',domnode));
    % Sensor specific portions of metadata
    output_meta{i}.CollectionInfo.Parameter{1}.name = 'SLICE';
    output_meta{i}.CollectionInfo.Parameter{1}.value = num2str(slice);
    output_meta{i}.CollectionInfo.Parameter{2}.name = 'SWATH';
    output_meta{i}.CollectionInfo.Parameter{2}.value = swath;
    output_meta{i}.CollectionInfo.Parameter{3}.name = 'BURST';
    output_meta{i}.CollectionInfo.Parameter{3}.value = num2str(i);
    
    %% ImageData.ValidData
    % Get valid bounds of burst from metadata.  Assume a rectangular valid
    % area-- not totally true, but all that seems to be defined by the
    % product XML metadata.
    if num_bursts>0 % Valid data does not seem to be defined for STRIPMAP data...
        firstSamples = str2num(char(xp.evaluate(...
            ['product/swathTiming/burstList/burst[' num2str(i) ']/firstValidSample'],domnode)));
        lastSamples = str2num(char(xp.evaluate(...
            ['product/swathTiming/burstList/burst[' num2str(i) ']/lastValidSample'],domnode)));
        valid_cols = (firstSamples>=0)&(lastSamples>=0);
        first_col = find(valid_cols,1,'first') - 1; % SICD zero-based vs MATLAB one-based
        last_col = find(valid_cols,1,'last') - 1; % SICD zero-based vs MATLAB one-based
        first_row = max(firstSamples(valid_cols)); % Sentinel XML and SICD both zero-based, no -1 needed
        last_row = min(lastSamples(valid_cols)); % Sentinel XML and SICD both zero-based
        % From SICD spec: Vertices ordered clockwise with vertex 1
        % determined by: (1) minimum row index, (2) minimum column index if
        % 2 vertices exist with minimum row index.
        output_meta{i}.ImageData.ValidData.Vertex(1).Row = first_row;
        output_meta{i}.ImageData.ValidData.Vertex(1).Col = first_col;
        output_meta{i}.ImageData.ValidData.Vertex(2).Row = first_row;
        output_meta{i}.ImageData.ValidData.Vertex(2).Col = last_col;
        output_meta{i}.ImageData.ValidData.Vertex(3).Row = last_row;
        output_meta{i}.ImageData.ValidData.Vertex(3).Col = last_col;
        output_meta{i}.ImageData.ValidData.Vertex(4).Row = last_row;
        output_meta{i}.ImageData.ValidData.Vertex(4).Col = first_col;
    end
    
    %% Timeline
    if num_bursts>0
        % This is the first and last zero doppler times of the columns in
        % the burst.  This isn't really what we mean by CollectStart and
        % CollectDuration in SICD (really we want first and last pulse
        % times), but its all we have.
        [start_s, start_frac] = datenum_w_frac(char(xp.evaluate(...
            ['product/swathTiming/burstList/burst[' num2str(i) ']/azimuthTime'],domnode)));
        first_line_relative_start = 0; % CollectStart is zero Doppler time of first column
    else % STRIPMAP
        [start_s, start_frac] = datenum_w_frac(char(xp.evaluate(...
            'product/generalAnnotation/downlinkInformationList/downlinkInformation/firstLineSensingTime',domnode)));
        % Maybe CollectStart/CollectDuration should be set by
        % product/imageAnnotation/imageInformation/productFirstLineUtcTime
        % and productLastLineUtcTime.  This would make it consistent with
        % non-stripmap which just defines first and last zero doppler
        % times, but is not really consistent with what SICD generally
        % means by CollectStart/CollectDuration.
        [stop_s, stop_frac] = datenum_w_frac(char(xp.evaluate(...
            'product/generalAnnotation/downlinkInformationList/downlinkInformation/lastLineSensingTime',domnode)));
        output_meta{i}.Timeline.CollectStart = start_s + (start_frac/SECONDS_IN_A_DAY);
        output_meta{i}.Timeline.CollectDuration = ...
            round((stop_s - start_s)*SECONDS_IN_A_DAY) + (stop_frac-start_frac);
        first_line_relative_start = ... % Convert from UTC to relative to start
            round((azimuth_time_first_line-start_s)*SECONDS_IN_A_DAY) + ... % Convert to seconds
            (azimuth_time_frac-start_frac); % Handle fractional seconds
    end

    % After we have start_s, we can generate CoreName
    output_meta{i}.CollectionInfo.CoreName = [...
        ... % Prefix with the NGA CoreName standard format
        upper(datestr(start_s,'ddmmmyy')) ...
        common_meta.CollectionInfo.CollectorName ...
        ... % The following core name is unique within all Sentinel-1
        ...% coherent data periods:
        num2str(str2double(xp.evaluate(...
        'product/adsHeader/missionDataTakeId',domnode)),'%07u') '_' ... % Is 7 digits enough for the life of the mission?
        num2str(slice, '%02u') '_' ... % Slice
        swath '_' ... % Swath
        num2str(i, '%02u')]; % Burst
    
    %% Position
    % Polynomial is computed with respect to time from start of burst
    state_vector_T_burst = round((state_vector_T-start_s)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (state_vector_T_frac-start_frac); % Handle fractional seconds
    % sv2poly.m shows ways to determine best polynomial order, but 5th is almost always best
    polyorder = min(5, numel(state_vector_T_burst) - 1);
    P_x = polyfit(state_vector_T_burst, state_vector_pos(:,1), polyorder);
    P_y = polyfit(state_vector_T_burst, state_vector_pos(:,2), polyorder);
    P_z = polyfit(state_vector_T_burst, state_vector_pos(:,3), polyorder);
    output_meta{i}.Position.ARPPoly.X = P_x(end:-1:1).';
    output_meta{i}.Position.ARPPoly.Y = P_y(end:-1:1).';
    output_meta{i}.Position.ARPPoly.Z = P_z(end:-1:1).';

    %% RMA
    % Sentinel-1 is always right-looking, so TimeCAPoly should never have
    % to be "flipped" for left-looking cases.
    output_meta{i}.RMA.INCA.TimeCAPoly = first_line_relative_start + eta_mid; % SCP zero Doppler time relative to start
    output_meta{i}.RMA.INCA.TimeCAPoly(2,1) = ss_zd_s/common_meta.Grid.Col.SS; % Convert zero doppler spacing from sec/pixels to sec/meters

    %% Doppler centroid
    % We choose the single Doppler centroid polynomial closest to the
    % center of the current burst.
    dc_est_times = round((dc_az_time_s-start_s) * SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (dc_az_time_frac-start_frac); % Handle fractional seconds
    [~,dc_poly_ind] = min(abs(dc_est_times - output_meta{i}.RMA.INCA.TimeCAPoly(1)));
    % Shift polynomial from origin at dc_t0 (reference time for Sentinel
    % polynomial) to SCP time (reference time for SICD polynomial)
    range_time_scp = common_meta.RMA.INCA.R_CA_SCP*2/SPEED_OF_LIGHT;
    % The Doppler centroid field in the Sentinel-1 metadata is not
    % complete, so we cannot use it directly.  That description of Doppler
    % centroid by itself does not vary by azimuth although the
    % Col.DeltaKCOAPoly we see in the data definitely does. We will define
    % DopCentroidPoly differently later down in the code.  We also comment
    % out any definition of TimeCOAPoly based off of this DopCentroidPoly
    % definition, since it is wrong without some further adjustment.
    % dop_rate_poly_rg_shifted = polyshift(data_dc_poly{dc_poly_ind}(:), ...
    %     range_time_scp - dc_t0(dc_poly_ind));
    % % Scale 1D polynomial to from Hz/s^n to Hz/m^n
    % output_meta{i}.RMA.INCA.DopCentroidPoly = dop_rate_poly_rg_shifted.*...
    %     (2/SPEED_OF_LIGHT).^(0:(length(dop_rate_poly_rg_shifted)-1)).';

    %% Doppler rate
    % Total Doppler rate is a combination of the Doppler FM rate and the
    % Doppler rate introduced by the scanning of the antenna.
    % We pick a single velocity magnitude at closest approach to represent
    % the entire burst.  This is valid, since the magnitude of the velocity
    % changes very little.
    pos_coefs = [P_x(:) P_y(:) P_z(:)];
    % Velocity is derivate of position.
    vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
    vel_x = polyval(vel_coefs(:,1), output_meta{i}.RMA.INCA.TimeCAPoly(1));
    vel_y = polyval(vel_coefs(:,2), output_meta{i}.RMA.INCA.TimeCAPoly(1));
    vel_z = polyval(vel_coefs(:,3), output_meta{i}.RMA.INCA.TimeCAPoly(1));
    vm_ca = sqrt(vel_x.^2 + vel_y.^2 + vel_z.^2); % Magnitude of the velocity at SCP closest approach
    % Compute FM Doppler Rate, k_a
    % We choose the single azimuth FM rate polynomial closest to the
    % center of the current burst.
    az_rate_times = round((az_t_s-start_s) * SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (az_t_frac-start_frac); % Handle fractional seconds
    [~,az_rate_poly_ind] = min(abs(az_rate_times - output_meta{i}.RMA.INCA.TimeCAPoly(1)));
    % SICD's Doppler rate seems to be FM Doppler rate, not total Doppler rate
    % Shift polynomial from origin at az_t0 (reference time for Sentinel
    % polynomial) to SCP time (reference time for SICD polynomial)
    DR_CA = polyshift(k_a_poly{az_rate_poly_ind}(:), ...
        range_time_scp - az_t0(az_rate_poly_ind));
    % Scale 1D polynomial to from Hz/s^n to Hz/m^n
    DR_CA = DR_CA.*(2/SPEED_OF_LIGHT).^(0:(numel(DR_CA)-1)).';
    % Polynomial representing range as a function of range distance from SCP
    r_ca = [output_meta{i}.RMA.INCA.R_CA_SCP; 1];
    % RMA.INCA.DRateSFPoly is a function of Doppler rate.
    output_meta{i}.RMA.INCA.DRateSFPoly = - conv(DR_CA, r_ca) * ...
        SPEED_OF_LIGHT / (2 * fc * (vm_ca^2)); % Assumes a SGN of -1
    
    %% TimeCOAPoly
    % TimeCOAPoly = TimeCA + (DopCentroid/dop_rate); % True if DopCentroidCOA = true
    % Since we don't know how to evaluate this equation analytically, we
    % could evaluate samples of it across our image and fit a 2D polynomial
    % to it later.
    POLY_ORDER = 2; % Same order of native metadata polynomials to describe Doppler Centroid, azimuth fm rate, etc.
    grid_samples = POLY_ORDER + 1; % in each dimension
    % For debugging, this lets us compute phase for every pixel:
    % cols = linspace(1, double(output_meta{i}.ImageData.NumCols), output_meta{i}.ImageData.NumCols);
    % rows = linspace(1, double(output_meta{i}.ImageData.NumRows), output_meta{i}.ImageData.NumRows);
    cols = round(linspace(1, double(output_meta{i}.ImageData.NumCols), grid_samples));
    rows = round(linspace(1, double(output_meta{i}.ImageData.NumRows), grid_samples));
    coords_az_m = double(cols - (output_meta{i}.ImageData.SCPPixel.Col+1)) * ...
        output_meta{i}.Grid.Col.SS;
    coords_rg_m = double(rows - (output_meta{i}.ImageData.SCPPixel.Row+1)) * ...
        output_meta{i}.Grid.Row.SS;
    timeca_sampled = sicd_polyval2d(output_meta{i}.RMA.INCA.TimeCAPoly(:).',coords_az_m,coords_rg_m);
    % dopcentroid_sampled = sicd_polyval2d(output_meta{i}.RMA.INCA.DopCentroidPoly,coords_az_m,coords_rg_m);
    doprate_sampled = sicd_polyval2d(DR_CA,coords_az_m,coords_rg_m);
    % timecoa_sampled = timeca_sampled + (dopcentroid_sampled./doprate_sampled);
    
    %% Grid.Col.DeltaKCOAPoly
    % Reference: Definition of the TOPS SLC deramping function for products
    % generated by the S-1 IPF, COPE-GSEG-EOPG-TN-14-0025
    tau = tau_0 + delta_tau_s * (0:double(common_meta.ImageData.NumRows-1)); % Range time for each sample
    % The vm_ca used here is slightly different than the ESA deramp
    % document, since the document interpolates the velocity values given
    % rather than the position values, which is what we do here.
    k_s = (2 * vm_ca / SPEED_OF_LIGHT) * fc * k_psi; % Doppler rate as introduced by antenna steering, k_s
    k_a = polyval(k_a_poly{az_rate_poly_ind}(end:-1:1), tau - az_t0(az_rate_poly_ind)); % Doppler FM rate
    k_t = (k_a * k_s)./(k_a - k_s); % Total Doppler Centroid Rate
    f_eta_c= polyval(data_dc_poly{dc_poly_ind}(end:-1:1), tau - dc_t0(dc_poly_ind)); % Doppler Centroid
    % This old implementation was imprecise.  The limits were one more than
    % the spacing.  Should have been a -1 in here somewhere.  New
    % formulation is clearer (at least to me). -WCS
    % eta = linspace(-double(output_meta{i}.ImageData.NumCols)*ss_zd_s/2, ...
    %     double(output_meta{i}.ImageData.NumCols)*ss_zd_s/2, ...
    %     double(output_meta{i}.ImageData.NumCols));
    eta = (-double(output_meta{i}.ImageData.SCPPixel.Col) * ss_zd_s) + ...
        (0:double(output_meta{i}.ImageData.NumCols-1))*ss_zd_s;
    eta_c = - f_eta_c./k_a; % Beam center crossing time.  TimeCOA in SICD terminology
    eta_ref = eta_c - eta_c(1);
    [eta_grid, eta_ref_grid] = ndgrid(eta(cols), eta_ref(rows));
    eta_arg = eta_grid - eta_ref_grid;
    k_t_grid = repmat(k_t(rows),numel(cols),1);
    f_eta_c_grid = repmat(f_eta_c(rows),numel(cols),1);
    deramp_phase = k_t_grid.*(eta_arg.^2)/2;
    demod_phase = f_eta_c_grid.*eta_arg;
    total_phase = deramp_phase + demod_phase; % Sampled phase correction for deramping and demodding

    %% Least squares fit for 2D polynomials
    % A*x = b
    [coords_az_m_2d, coords_rg_m_2d] = ndgrid(coords_az_m, coords_rg_m);
    a = zeros(numel(total_phase), (POLY_ORDER+1)^2);
    for k = 0:POLY_ORDER
        for j = 0:POLY_ORDER
            a(:,k*(POLY_ORDER+1)+j+1) = (coords_rg_m_2d(:).^j).*(coords_az_m_2d(:).^k);
        end
    end
    % b_coa = zeros((POLY_ORDER+1)^2,1);
    b_phase = zeros((POLY_ORDER+1)^2,1);
    for k=1:((POLY_ORDER+1)^2)
       % b_coa(k)=sum(timecoa_sampled(:).*a(:,k)); % center of aperture
       b_phase(k)=sum(total_phase(:).*a(:,k)); % phase to deramp in azimuth
    end
    A=zeros((POLY_ORDER+1)^2);
    for k=1:((POLY_ORDER+1)^2)
        for j=1:((POLY_ORDER+1)^2)
            A(k,j)=sum(a(:,k).*a(:,j));
        end
    end
    % MATLAB often flags these as badly scaled, but results still appear valid
    old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
    % x_coa=A\b_coa;
    x_phase=A\b_phase;
    warning(old_warning_state);
    phase=reshape(x_phase, POLY_ORDER+1, POLY_ORDER+1);
    % output_meta{i}.Grid.TimeCOAPoly=reshape(x_coa, POLY_ORDER+1, POLY_ORDER+1);
    % DeltaKCOAPoly is derivative of phase in Col direction
    output_meta{i}.Grid.Col.DeltaKCOAPoly = phase(:,2:end) .* ...
        repmat((1:(size(phase,2)-1)),size(phase,1),1);
    
    %% DopCentroidPoly/TimeCOAPoly
    % Another way to derive the Doppler Centroid, which is back-calculated
    % from the ESA-documented azimuth deramp phase function.
    output_meta{i}.RMA.INCA.DopCentroidPoly = output_meta{i}.Grid.Col.DeltaKCOAPoly * ...
        common_meta.Grid.Col.SS / ss_zd_s;
    dopcentroid2_sampled = sicd_polyval2d(output_meta{i}.RMA.INCA.DopCentroidPoly, coords_az_m, coords_rg_m);
    timecoa_sampled = timeca_sampled + (dopcentroid2_sampled./doprate_sampled);
    % Convert sampled TimeCOA to polynomial
    b_coa = zeros((POLY_ORDER+1)^2,1);
    for k=1:((POLY_ORDER+1)^2)
       b_coa(k)=sum(timecoa_sampled(:).*a(:,k)); % center of aperture
    end
    old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
    x_coa=A\b_coa;
    warning(old_warning_state);
    output_meta{i}.Grid.TimeCOAPoly=reshape(x_coa, POLY_ORDER+1, POLY_ORDER+1);

    %% Timeline
    % We don't know the precise start and stop time of each burst (as in
    % the times of first and last pulses), so we use the min and max COA
    % time, which is a closer approximation than the min and max zero
    % Doppler times.  At least COA times will not overlap between bursts.
    if num_bursts>0 % STRIPMAP case uses another time origin different from first zero Doppler time
        time_offset = min(timecoa_sampled(:));
        % Datetime more precise than serial date number.  This precision is
        % important for TOPSAR burst relationships.  Datetime introduced in
        % MATLAB 2014b.
        output_meta{i}.Timeline.CollectStart = ... % Without fractional seconds
            datetime(start_s,'ConvertFrom','datenum');
        output_meta{i}.Timeline.CollectStart.Second = ... % With fraction seconds
            output_meta{i}.Timeline.CollectStart.Second + start_frac + time_offset;
        output_meta{i}.Timeline.CollectDuration = ...
            max(timecoa_sampled(:)) - min(timecoa_sampled(:));
        % Adjust all SICD fields that were dependent on start time
        % Time is output of polynomial:
        output_meta{i}.Grid.TimeCOAPoly(1) = ...
            output_meta{i}.Grid.TimeCOAPoly(1) - time_offset;
        output_meta{i}.RMA.INCA.TimeCAPoly(1) = ...
            output_meta{i}.RMA.INCA.TimeCAPoly(1) - time_offset;
        % Time is input of polynomial:
        for j = 'XYZ'
            output_meta{i}.Position.ARPPoly.(j) = ...
                polyshift(output_meta{i}.Position.ARPPoly.(j),time_offset);
        end
    end
    output_meta{i}.Timeline.IPP.Set.TStart=0;
    output_meta{i}.Timeline.IPP.Set.TEnd=output_meta{i}.Timeline.CollectDuration;
    output_meta{i}.Timeline.IPP.Set.IPPStart=uint32(0);
    output_meta{i}.Timeline.IPP.Set.IPPEnd=uint32(floor(output_meta{i}.Timeline.CollectDuration*prf)); % An approximation
    output_meta{i}.ImageFormation.TEndProc=output_meta{i}.Timeline.CollectDuration;

    %% GeoData
    % Rough estimate of SCP (interpolated from metadata geolocation grid)
    % to bootstrap (point_image_to_ground uses it only to find tangent to
    % ellipsoid.)  Then we will immediately replace it with a more precise
    % value from point_image_to_ground and the SICD sensor model.
    output_meta{i}.GeoData.SCP.LLH.Lat = scp_lats(i);
    output_meta{i}.GeoData.SCP.LLH.Lon = scp_lons(i);
    output_meta{i}.GeoData.SCP.LLH.HAE = scp_hgts(i);
    % Note that blindly using the heights in the geolocationGridPointList
    % can result in some confusing results.  Since the scenes can be
    % extremely large, you could easily be using a height in your
    % geolocationGridPointList that is very high, but still have ocean
    % shoreline in your scene. Blindly projecting the the plane tangent to
    % the inflated ellipsoid at SCP could result in some badly placed
    % geocoords in Google Earth.  Of course, one must always be careful
    % with ground projection and height variability, but probably even more
    % care is warranted in this data than even usual due to large scene
    % sizes and frequently steep graze angles.
    % Note also that some Sentinel-1 data we have see has different heights
    % in the geolocation grid for polarimetric channels from the same
    % swath/burst!?!
    ecf=geodetic_to_ecf([output_meta{i}.GeoData.SCP.LLH.Lat ...
        output_meta{i}.GeoData.SCP.LLH.Lon output_meta{i}.GeoData.SCP.LLH.HAE]);
    output_meta{i}.GeoData.SCP.ECF.X=ecf(1);
    output_meta{i}.GeoData.SCP.ECF.Y=ecf(2);
    output_meta{i}.GeoData.SCP.ECF.Z=ecf(3);
    % Now that SCP has been populated, populate GeoData.SCP more precisely.
    ecf = point_image_to_ground([common_meta.ImageData.SCPPixel.Row; common_meta.ImageData.SCPPixel.Col], output_meta{i});
    output_meta{i}.GeoData.SCP.ECF.X=ecf(1);
    output_meta{i}.GeoData.SCP.ECF.Y=ecf(2);
    output_meta{i}.GeoData.SCP.ECF.Z=ecf(3);
    llh=ecf_to_geodetic([output_meta{i}.GeoData.SCP.ECF.X output_meta{i}.GeoData.SCP.ECF.Y output_meta{i}.GeoData.SCP.ECF.Z]);
    output_meta{i}.GeoData.SCP.LLH.Lat=llh(1);
    output_meta{i}.GeoData.SCP.LLH.Lon=llh(2);
    output_meta{i}.GeoData.SCP.LLH.HAE=llh(3);

    %% SCPCOA (and other stuff)
    output_meta{i} = derived_sicd_fields(output_meta{i});
end

end


% Reproduces hamming functionality from the MATLAB Signal Processing
% Toolbox, but allows for arbitrary coefficients of raised cosine.
function w = raised_cos_fun(n, coef)
    w = coef - (1-coef)*cos(2*pi*(0:ceil(n/2)-1)'/(n-1));
    if ~rem(n,2)
        w = [w; w(end:-1:1)];
    else
        w = [w; w(end-1:-1:1)];
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

% Works on both restituted orbit files and precise orbit ephemerides orbit
% files.  Time range is a two-element vector of MATLAB date numbers.
% Returned state vectors will be within a second of these bounds.
function [ times_sec, times_frac, pos, vel ] = get_osv_from_eof( eof_domnode, time_range )

SECONDS_IN_A_DAY = 24*60*60;

osv_list = eof_domnode.getElementsByTagName('OSV');
[times_sec, times_frac] = deal(zeros(osv_list.getLength,1));
[pos, vel] = deal(zeros(osv_list.getLength,3));
valid = false(osv_list.getLength,1);
for i = 1:osv_list.getLength
    % Orbit files have 3 times (TAI, UTC, and UT1)
    % We choose UTC here, since 1) that is what is used in the orbit state
    % vectors in the SLC annotation file and 2) that is the SICD standard.
    % Restituted orbit files seem to have state vectors as the exact same
    % UTC times as those in the state vectors in the SLC annotation files,
    % while precise orbit files have state vectors on whole UTC seconds.
    time_reference = 'UTC';
    time_date_str = char(osv_list.item(i-1).getElementsByTagName(time_reference).item(0).getFirstChild().getData());
    time_date_str(1:strfind(time_date_str,'='))=''; % Remove 'UTC='
    [times_sec(i), times_frac(i)]= datenum_w_frac(time_date_str);
    if (nargin < 2) || ... % If no time ranges were passed in, read in all
            (((times_sec(i) - time_range(2)) - (1/SECONDS_IN_A_DAY) <= 0) && ...
            ((times_sec(i) - time_range(1)) + (1/SECONDS_IN_A_DAY) >= 0))
        pos(i,1) = str2double(osv_list.item(i-1).getElementsByTagName('X').item(0).getFirstChild().getData());
        pos(i,2) = str2double(osv_list.item(i-1).getElementsByTagName('Y').item(0).getFirstChild().getData());
        pos(i,3) = str2double(osv_list.item(i-1).getElementsByTagName('Z').item(0).getFirstChild().getData());
        vel(i,1) = str2double(osv_list.item(i-1).getElementsByTagName('VX').item(0).getFirstChild().getData());
        vel(i,2) = str2double(osv_list.item(i-1).getElementsByTagName('VY').item(0).getFirstChild().getData());
        vel(i,3) = str2double(osv_list.item(i-1).getElementsByTagName('VZ').item(0).getFirstChild().getData());
        valid(i) = true;
        if ~strcmpi('NOMINAL', osv_list.item(i-1).getElementsByTagName('Quality').item(0).getFirstChild().getData())
            warning('META2SICD_S1PRODUCT:DEGRADED_EPHEMERIS', ...
                'Ephemeris found in requested file is of degraded quality.');
        end
    end
end
if ~any(valid)
    error('META2SICD_S1PRODUCT:NON_OVERLAPPING_EOF','EOF file does not span the same times as SLC.');
end
pos = pos(valid,:);
vel = vel(valid,:);
times_sec = times_sec(valid);
times_frac = times_frac(valid);

% XPath version.  Goes MUCH slower!  Don't use.
% xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
% base_xpath = 'Earth_Explorer_File/Data_Block/List_of_OSVs/OSV';
% num_grid_points=str2double(xp.evaluate(['count(' base_xpath ')'],dom_node));
% for i = 1:num_grid_points
%     time_date_str = char(xp.evaluate([base_xpath '[' num2str(i) ']/UTC'],dom_node));
%     time_date_str(1:strfind(time_date_str,'='))=''; % Remove 'UTC='
%     [times_sec(i), times_frac(i)]= datenum_w_frac(time_date_str);
%     pos(i,1) = str2double(xp.evaluate([base_xpath '[' num2str(i) ']/X'],dom_node));
%     pos(i,2) = str2double(xp.evaluate([base_xpath '[' num2str(i) ']/Y'],dom_node));
%     pos(i,3) = str2double(xp.evaluate([base_xpath '[' num2str(i) ']/Z'],dom_node));
% end

end


% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////