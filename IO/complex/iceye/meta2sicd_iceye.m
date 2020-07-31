function [ output_meta ] = meta2sicd_iceye(hdffile)
%META2SICD_ICEYE Converts ICEYE HDF5 into a SICD metadata structure
%
% Written by: Tim Cox, NRL; Wade Schwartzkopf, NGA/Reseach
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

SECONDS_IN_A_DAY = 24*60*60;

%% Collection Info
output_meta.CollectionInfo.CollectorName = get_hdf_data(hdffile,'/','satellite_name'); 
output_meta.CollectionInfo.CoreName = get_hdf_data(hdffile,'/','product_name');
output_meta.CollectionInfo.CollectType='MONOSTATIC';
output_meta.CollectionInfo.RadarMode.ModeType = upper(get_hdf_data(hdffile,'/','acquisition_mode'));
output_meta.CollectionInfo.RadarMode.ModeID = get_hdf_data(hdffile,'/','product_type');
output_meta.CollectionInfo.Classification = 'UNCLASSIFIED';

%% Image Creation
output_meta.ImageCreation.Application = ['ICEYE_P_' ...
    num2str(get_hdf_data(hdffile,'/','processor_version'))]; 
output_meta.ImageCreation.DateTime = datenum(reshape(get_hdf_data(...
    hdffile,'/','processing_time'),1,[]),'yyyy-mm-ddTHH:MM:SS'); 
output_meta.ImageCreation.Profile='Prototype';

%% Image Data
switch get_hdf_data(hdffile,'/','sample_precision')
    case 'float32'
        output_meta.ImageData.PixelType='RE32F_IM32F';
    case 'int16'
        output_meta.ImageData.PixelType='RE16I_IM16I';
    otherwise
        warning('META2SICD_ICEYE:UnrecognizedSamplePrecision', 'Only float32 and int16 precisions supported.');
end
output_meta.ImageData.NumRows = uint32(get_hdf_data(hdffile,'/','number_of_range_samples')); 
output_meta.ImageData.NumCols = uint32(get_hdf_data(hdffile,'/','number_of_azimuth_samples'));
output_meta.ImageData.FirstRow = 0;
output_meta.ImageData.FirstCol = 0;
output_meta.ImageData.FullImage.NumRows = output_meta.ImageData.NumRows;
output_meta.ImageData.FullImage.NumCols = output_meta.ImageData.NumCols;
coord_center = get_hdf_data(hdffile,'/','coord_center');
output_meta.ImageData.SCPPixel.Row = uint32(coord_center(1))-1;
if coord_center(2)>0 && coord_center(2)<output_meta.ImageData.NumCols
    output_meta.ImageData.SCPPixel.Col = uint32(coord_center(2))-1;
    if strcmpi(get_hdf_data(hdffile,'/','look_side'),'left')
        output_meta.ImageData.SCPPixel.Col = output_meta.ImageData.NumCols - ...
            output_meta.ImageData.SCPPixel.Col - 1;
    end
else
    % Bug in earlier version of ICEYE processor sometimes resulted in
    % center pixel outside of image in the column direction.  In this case,
    % we just choose our own.  The choice of SCP is essentially arbitrary
    % anyway.
    output_meta.ImageData.SCPPixel.Col = output_meta.ImageData.NumCols/2;
end

%% GeoData
avg_scene_height = double(get_hdf_data(hdffile,'/','avg_scene_height')); 
SCP = geodetic_to_ecf(coord_center(3),coord_center(4),avg_scene_height);
if strcmpi(get_hdf_data(hdffile,'/','geo_ref_system'),'WGS84')
    output_meta.GeoData.EarthModel = 'WGS_84';  % Should always be the case
else
    warning('META2SICD_ICEYE:UnexpectedEarthModel', 'Only WGS84 expected for earth model.');
end
output_meta.GeoData.SCP.ECF.X = SCP(1);
output_meta.GeoData.SCP.ECF.Y = SCP(2);
output_meta.GeoData.SCP.ECF.Z = SCP(3);
output_meta.GeoData.SCP.LLH.Lat = coord_center(3);
output_meta.GeoData.SCP.LLH.Lon = coord_center(4);
output_meta.GeoData.SCP.LLH.HAE = avg_scene_height;
% Image corner coordinates
% We could do it this way, but we would have to swap first/last for left
% looking cases.  Instead, we will let derived_sicd_fields compute these.
% FRFC = get_hdf_data(hdffile,'/','coord_first_near'); 
% output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=FRFC(3);
% output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=FRFC(4);
% FRLC = get_hdf_data(hdffile,'/','coord_last_near'); 
% output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=FRLC(3);
% output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=FRLC(4);
% LRFC = get_hdf_data(hdffile,'/','coord_first_far'); 
% output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=LRFC(3);
% output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=LRFC(4);
% LRLC = get_hdf_data(hdffile,'/','coord_last_far'); 
% output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=LRLC(3);
% output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=LRLC(4);


%% Grid
% We will derive Grid.Col.SS/ImpRespBW later
% These two computations of Grid.Row.SS should be identical
% output_meta.Grid.Row.SS = get_hdf_data(hdffile,'/','slant_range_spacing');
output_meta.Grid.Row.SS = SPEED_OF_LIGHT/(2*get_hdf_data(hdffile,'/','range_sampling_rate'));
fc = get_hdf_data(hdffile,'/','carrier_frequency'); %Hz
output_meta.Grid.Row.KCtr=2*fc/SPEED_OF_LIGHT;
output_meta.Grid.Col.KCtr=0;
output_meta.Grid.Type='RGZERO'; %based on ICEYE spec, all SLC files are range, zero Doppler
output_meta.Grid.ImagePlane = 'SLANT'; %SLC is always slant
output_meta.Grid.Row.Sgn=-1; 
output_meta.Grid.Col.Sgn=-1;
% ICEYE documentation says chirp_bandwidth should be float, but it is int
% in datasets we have seen.
output_meta.Grid.Row.ImpRespBW = 2*double(get_hdf_data(hdffile,'/','chirp_bandwidth'))/SPEED_OF_LIGHT;
output_meta.Grid.Row.DeltaKCOAPoly=0;
output_meta.Grid.Row.WgtType.WindowName = get_hdf_data(hdffile,'/','window_function_range');
output_meta.Grid.Col.WgtType.WindowName = get_hdf_data(hdffile,'/','window_function_azimuth');
% For now only have seen uniform weighted data
if strcmpi(output_meta.Grid.Row.WgtType.WindowName,'NONE')
    output_meta.Grid.Row.WgtType.WindowName = 'UNIFORM';
end
% Although no weighting is applied, ICEYE generally does not remove
% slow-time beam shape, so spectrum will not necessarily be flat.  (This is
% reflected in the STBeamComp parameter in SICD.)
if strcmpi(output_meta.Grid.Col.WgtType.WindowName,'NONE')
    output_meta.Grid.Col.WgtType.WindowName = 'UNIFORM';
end
% ImpRespWid will be computed in derived_sicd_fields later

%% Timeline
[start_t, start_frac] = datenum_w_frac(reshape(get_hdf_data(hdffile,'/','acquisition_start_utc'),1,[]));
[end_t, end_frac] = datenum_w_frac(reshape(get_hdf_data(hdffile,'/','acquisition_end_utc'),1,[]));
output_meta.Timeline.CollectStart = start_t + (start_frac/SECONDS_IN_A_DAY);
output_meta.Timeline.CollectDuration = round((end_t-start_t)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (end_frac-start_frac); % Handle fractional seconds
output_meta.Timeline.IPP.Set.TStart = 0; 
output_meta.Timeline.IPP.Set.TEnd = output_meta.Timeline.CollectDuration;
acquisition_prf = get_hdf_data(hdffile,'/','acquisition_prf');
output_meta.Timeline.IPP.Set.IPPStart = 0;
% We don't realy know how many IPPs.  Just making data consistent
output_meta.Timeline.IPP.Set.IPPEnd = round(acquisition_prf*...
    (output_meta.Timeline.IPP.Set.TEnd-output_meta.Timeline.IPP.Set.TStart));
output_meta.Timeline.IPP.Set.IPPPoly = [0; acquisition_prf]; %assume constant PRF

%% Position
state_vector_time_str = get_hdf_data(hdffile,'/','state_vector_time_utc').'; 
state_vector_T  = zeros(1,size(state_vector_time_str,1));
state_vector_T_frac  = zeros(1,size(state_vector_time_str,1));
for i=1:size(state_vector_time_str,1)
    [state_vector_T(i), state_vector_T_frac(i)] = datenum_w_frac(state_vector_time_str(i,:));
end
% Need to make this polynomial relative to the start time of the collect
t = round((state_vector_T-start_t)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
        (state_vector_T_frac-start_frac); % Handle fractional seconds
output_meta.Position.ARPPoly.X = fliplr(polyfit(t(:),get_hdf_data(hdffile,'/','posX'),5))';
output_meta.Position.ARPPoly.Y = fliplr(polyfit(t(:),get_hdf_data(hdffile,'/','posY'),5))';
output_meta.Position.ARPPoly.Z = fliplr(polyfit(t(:),get_hdf_data(hdffile,'/','posZ'),5))';
% Velocity is in native metadata here, but we will derive it
% velX = get_hdf_data(hdffile,'/','velX');
% velY = get_hdf_data(hdffile,'/','velY');
% velZ = get_hdf_data(hdffile,'/','velZ');

%% RadarCollection
Pol = get_hdf_data(hdffile,'/','polarization');
output_meta.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization = [Pol(1) ':' Pol(2)];
output_meta.RadarCollection.TxPolarization = Pol(1);
TxBandwidth = double(get_hdf_data(hdffile,'/','chirp_bandwidth'));  % Hz
output_meta.RadarCollection.TxFrequency.Min = fc-TxBandwidth/2;
output_meta.RadarCollection.TxFrequency.Max = fc+TxBandwidth/2;
output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart = ...
    output_meta.RadarCollection.TxFrequency.Min;
output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth = TxBandwidth;
output_meta.RadarCollection.Waveform.WFParameters.TxPulseLength = ...
    get_hdf_data(hdffile,'/','chirp_duration');  % sec
output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate = ...
    get_hdf_data(hdffile,'/','range_sampling_rate');  % Hz
output_meta.RadarCollection.Waveform.WFParameters.RcvDemodType = 'CHIRP';
output_meta.RadarCollection.Waveform.WFParameters.RcvFMRate = 0;

%% ImageFormation
output_meta.ImageFormation.RcvChanProc=struct('NumChanProc',1,'PRFScaleFactor',1,'ChanIndex',1);
output_meta.ImageFormation.TxRcvPolarizationProc=[Pol(1) ':' Pol(2)];
output_meta.ImageFormation.ImageFormAlgo='RMA';
output_meta.ImageFormation.TStartProc=0;
output_meta.ImageFormation.TEndProc = output_meta.Timeline.CollectDuration;
output_meta.ImageFormation.TxFrequencyProc.MinProc = ...
    output_meta.RadarCollection.TxFrequency.Min;
output_meta.ImageFormation.TxFrequencyProc.MaxProc = ...
    output_meta.RadarCollection.TxFrequency.Max;
output_meta.ImageFormation.STBeamComp = 'NO';
output_meta.ImageFormation.ImageBeamComp = 'SV';
output_meta.ImageFormation.AzAutofocus='NO';
output_meta.ImageFormation.RgAutofocus='NO';

%% Radiometric
output_meta.Radiometric.BetaZeroSFPoly = get_hdf_data(hdffile,'/','calibration_factor');

%% RMA 
output_meta.RMA.RMAlgoType='OMEGA_K';
output_meta.RMA.ImageType='INCA';
% These two computations of near_range are the same:
% near_range = get_hdf_data(hdffile,'/','slant_range_to_first_pixel'); 
near_range = get_hdf_data(hdffile,'/','first_pixel_time')*SPEED_OF_LIGHT/2;
output_meta.RMA.INCA.R_CA_SCP = near_range + ...
    (double(output_meta.ImageData.SCPPixel.Row)*output_meta.Grid.Row.SS);
output_meta.RMA.INCA.FreqZero=fc;
ss_zd_s = get_hdf_data(hdffile,'/','azimuth_time_interval');
[zd_first, zd_first_frac]=datenum_w_frac(reshape(get_hdf_data(hdffile,'/','zerodoppler_start_utc'),1,[]));
[zd_last, zd_last_frac]=datenum_w_frac(reshape(get_hdf_data(hdffile,'/','zerodoppler_end_utc'),1,[]));
if strcmpi(get_hdf_data(hdffile,'/','look_side'),'left')
    ss_zd_s = -ss_zd_s;
    zd_left = zd_last;
    zd_left_frac = zd_last_frac;
else
    zd_left = zd_first;
    zd_left_frac = zd_first_frac;
end
% Zero doppler time of SCP relative to collect start
zd_t_scp = round((zd_left-start_t)*SECONDS_IN_A_DAY) + ... % Convert days to seconds
    (zd_left_frac-start_frac) + ... % Handle fractional seconds
    (double(output_meta.ImageData.SCPPixel.Col) * ss_zd_s);
% Reference time for Doppler polynomials according to ICEYE documentation
tref = get_hdf_data(hdffile,'/','first_pixel_time') + ...
    double(get_hdf_data(hdffile,'/','number_of_range_samples'))/...
    (2*get_hdf_data(hdffile,'/','range_sampling_rate'));
% Compute DRateSFPoly
pos_coefs = [output_meta.Position.ARPPoly.X(end:-1:1) ...
    output_meta.Position.ARPPoly.Y(end:-1:1) ...
    output_meta.Position.ARPPoly.Z(end:-1:1)];
% Velocity is derivate of position.
vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
vel_x = polyval(vel_coefs(:,1), zd_t_scp);
vel_y = polyval(vel_coefs(:,2), zd_t_scp);
vel_z = polyval(vel_coefs(:,3), zd_t_scp);
vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
r_ca = [output_meta.RMA.INCA.R_CA_SCP; 1]; % Polynomial representing range as a function of range distance from SCP
% Shift 1D polynomial to account for SCP (since we picked SCP to match
% reference point, this shouldn't do anything)
dop_rate_coeffs = get_hdf_data(hdffile,'/','doppler_rate_coeffs');
% Prior to ICEYE_P_1.14 processor, absolute value of Doppler rate was
% provided, not true Doppled rate.  Since Doppler rate is always negative
% will will flip sign appropriately.
dop_rate_coeffs = -sign(dop_rate_coeffs(1))* dop_rate_coeffs;
dop_rate_poly_rg_shifted=polyshift(dop_rate_coeffs, ...
    (2*output_meta.RMA.INCA.R_CA_SCP/SPEED_OF_LIGHT-tref));  % This offset should be zero because of how we have chosen SCP
% Scale 1D polynomial to from Hz/s^n to Hz/m^n
dop_rate_poly_rg_scaled=dop_rate_poly_rg_shifted.*...
    (2/SPEED_OF_LIGHT).^(0:(length(dop_rate_poly_rg_shifted)-1)).';
output_meta.RMA.INCA.DRateSFPoly = - conv(dop_rate_poly_rg_scaled,r_ca) * ... % Multiplication of two polynomials is just a convolution of their coefficients
    SPEED_OF_LIGHT / (2 * fc * vm_ca_sq(1)); % Assumes a SGN of -1
% Fields dependent on Doppler rate
% This computation of SS is actually better than the claimed SS
% (azimuth_ground_spacing) in many ways, because this makes all of the metadata
% internally consistent.  This must be the sample spacing exactly at SCP
% (which is the definition for SS in SICD), if the other metadata from
% which is it computed is correct and consistent. Since column SS can vary
% slightly over a RGZERO image, we don't know if the claimed sample spacing
% in the native metadata is at our chosen SCP, or another point, or an
% average across image or something else.
% output_meta.Grid.Col.SS = get_hdf_data(hdffile,'/','azimuth_ground_spacing');
output_meta.Grid.Col.SS = sqrt(vm_ca_sq(1)) * abs(ss_zd_s) * ...
    output_meta.RMA.INCA.DRateSFPoly(1,1);
% Convert to azimuth spatial bandwidth (cycles per meter)
dop_bw = get_hdf_data(hdffile,'/','total_processed_bandwidth_azimuth');
output_meta.Grid.Col.ImpRespBW = dop_bw*abs(ss_zd_s)/output_meta.Grid.Col.SS;
output_meta.RMA.INCA.TimeCAPoly = [zd_t_scp; ss_zd_s/output_meta.Grid.Col.SS];
% Doppler Centroid
dc_estimate_coeffs = get_hdf_data(hdffile,'/','dc_estimate_coeffs');
dc_estimate_time_utc = get_hdf_data(hdffile,'/','dc_estimate_time_utc').';
[dc_t, dc_t_frac, dc_zd_t]  = deal(zeros(1,size(dc_estimate_time_utc,1)));
for i=1:size(dc_estimate_time_utc,1)
    [dc_t(i), dc_t_frac(i)] = datenum_w_frac(dc_estimate_time_utc(i,:));
    dc_zd_t(i) = round((dc_t(i)-start_t)*SECONDS_IN_A_DAY) + ... % Convert days to seconds
        (dc_t_frac(i)-start_frac);  % Handle fractional seconds
end
M=49; % Arbitrary.  Seems as good as any.
f_dc=zeros(M,numel(dc_zd_t));
diff_t_rg = get_hdf_data(hdffile,'/','first_pixel_time') + linspace(0, ...
    double(get_hdf_data(hdffile,'/','number_of_range_samples'))/...
    get_hdf_data(hdffile,'/','range_sampling_rate'), M) - tref;
for n=1:numel(dc_zd_t)
    % Sampled values of Doppler centroid
    f_dc(:,n)=polyval(dc_estimate_coeffs(end:-1:1,n),diff_t_rg);
end
t_coa = repmat(dc_zd_t(:)', M, 1)+(f_dc./dop_rate_poly_rg_scaled(1,1));

% For SICD, we need polynomial as a function of range/azimuth distance from SCP
% in meters.
t_rg_scp = get_hdf_data(hdffile,'/','first_pixel_time') + ...
    double(output_meta.ImageData.SCPPixel.Row)/...
    get_hdf_data(hdffile,'/','range_sampling_rate');
range_scp_m = repmat((diff_t_rg(:) + tref - t_rg_scp) * (SPEED_OF_LIGHT/2), 1, numel(dc_zd_t));
azimuth_scp_m = repmat(output_meta.Grid.Col.SS * (dc_zd_t(:)' - zd_t_scp) / ss_zd_s, M, 1);

% Least squares fit (section 2.3 in paper)
% Compute the following continuous function for doppler centroid
% f_dc(azimuth_scp_m,range_scp_m)=...
%    x(1,1) + x(2,1)*range_scp_m + x(3,1)*range_scp_m*range_scp_m+...
%    x(1,2)*azimuth_scp_m + x(1,3)*azimuth_scp_m*azimuth_scp_m+...
%    x(2,2)*azimuth_scp_m*range_scp_m;
% A*x = b
% We add x(2,3),x(3,2), and x(3,3), even though they were not in the paper.
a=[ones(numel(f_dc),1) ...
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
output_meta.RMA.INCA.DopCentroidPoly=reshape(x,3,3);
output_meta.Grid.Col.DeltaKCOAPoly=...
    output_meta.RMA.INCA.DopCentroidPoly*ss_zd_s/output_meta.Grid.Col.SS;
output_meta.Grid.TimeCOAPoly=reshape(x2,3,3);
% Spotlight adjustments
if strcmp(output_meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    output_meta.Grid.TimeCOAPoly = output_meta.Grid.TimeCOAPoly(1);
    % This field required to compute Grid.Col.DeltaKCOAPoly, but not
    % allowed for spotlight in SICD.
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

%% SCPCOA
output_meta = derived_sicd_fields(output_meta);

end

% Reads HDF5 data
function value = get_hdf_data(hid_t,path,data_name)
gid = H5G.open(hid_t,path);
data = H5D.open(gid,data_name);
value = H5D.read(data);
% Versions of MATLAB prior to 2020a resulted in single-element cell arrays
% for some ICEYE HDF5 fields.
if iscell(value)
    value = value{1};
end
H5G.close(gid);
H5D.close(data);
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

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////