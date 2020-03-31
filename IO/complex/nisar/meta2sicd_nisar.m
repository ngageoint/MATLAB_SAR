function [ output_meta ] = meta2sicd_nisar( HDF5_fid )
%META2SICD_NISAR Converts NISAR SLC HDF5 into a SICD-style metadata structure
%
% Takes as input a file identifier to an open HDF5 file-- as is returned by
% H5F.open
%
% Note that MATLAB's h5disp is generally a good tool for manually browsing
% through HDF5 metadata like that found in NISAR format.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

SECONDS_IN_A_DAY = 24*60*60;

% Query for frequencies and polarization (each of which would create its
% own SICD)
freqs_id = H5D.open(HDF5_fid,'/science/LSAR/identification/listOfFrequencies');
freqs = H5D.read(freqs_id);
freqs = freqs(1,:);
H5D.close(freqs_id);
for i=1:numel(freqs)
    pol_id = H5D.open(HDF5_fid,['/science/LSAR/SLC/swaths/frequency' freqs(i) '/listOfPolarizations']);
    pols{i} = H5D.read(pol_id)';
    H5D.close(pol_id);
end

%% CollectionInfo
% CollectorName could also be deblank(get_hdf_data(HDF5_fid,'/science/LSAR/identification','missionId')')
output_meta.CollectionInfo.CollectorName=get_hdf_attribute(HDF5_fid,'/','mission_name');
% JPL suggested this as best way to form unique string for each collection
output_meta.CollectionInfo.CoreName=[...
    num2str(get_hdf_data(HDF5_fid,'/science/LSAR/identification','absoluteOrbitNumber'),'%07u') ...
    get_hdf_data(HDF5_fid,'/science/LSAR/identification','trackNumber')'];
output_meta.CollectionInfo.CollectType='MONOSTATIC';
output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
output_meta.CollectionInfo.Classification='UNCLASSIFIED';

%% ImageCreation
try
    % TODO: UNTESTED
    % Sample data did not have valid values for this
    output_meta.ImageCreation.Application=['ISCE ' ...
        deblank(get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/processingInformation/algorithms','ISCEVersion')')];
end
output_meta.ImageCreation.Profile='Prototype';

%% ImageData
output_meta.ImageData=struct(); % Just a placeholder
% Most subfields added below in "per band" section
% Used for computing SCP later

%% GeoData
output_meta.GeoData.EarthModel='WGS_84';
% Initially, we just seed this with a very rough value.  Later we will
% put in something more precise.
poly_str = get_hdf_data(HDF5_fid,'/science/LSAR/identification','boundingPolygon')';
poly_str(isletter(poly_str)|poly_str=='('|poly_str==')')='';
poly_str(poly_str==',')=';';
bound_poly=str2num(poly_str);
bound_poly=bound_poly(1:(end-1),:);
rough_center = mean(bound_poly);
output_meta.GeoData.SCP.LLH.Lat=rough_center(2);
output_meta.GeoData.SCP.LLH.Lon=rough_center(1);
output_meta.GeoData.SCP.LLH.HAE=...
    mean(get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/processingInformation/parameters/','referenceTerrainHeight'));
% NISAR stores this as single, which is sufficient for storage, but we need
% double precision for some of the computations later on, to include
% geopositioning of exact SCP.
output_meta.GeoData.SCP.LLH.HAE=double(output_meta.GeoData.SCP.LLH.HAE);
ecf=geodetic_to_ecf([output_meta.GeoData.SCP.LLH.Lat output_meta.GeoData.SCP.LLH.Lon output_meta.GeoData.SCP.LLH.HAE]);
output_meta.GeoData.SCP.ECF.X=ecf(1);
output_meta.GeoData.SCP.ECF.Y=ecf(2);
output_meta.GeoData.SCP.ECF.Z=ecf(3);
% Calling derived_sicd_fields at the end will populate these fields
% with the sensor model, so we don't need to do it here.
% latlon=get_hdf_attribute(dset_id(i),'Top Left Geodetic Coordinates');
% output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=bound_poly(1,1);
% output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=bound_poly(1,2);
% latlon=get_hdf_attribute(dset_id(i),'Bottom Left Geodetic Coordinates');
% output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=bound_poly(2,1);
% output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=bound_poly(2,2);
% latlon=get_hdf_attribute(dset_id(i),'Bottom Right Geodetic Coordinates');
% output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=bound_poly(3,1);
% output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=bound_poly(3,2);
% latlon=get_hdf_attribute(dset_id(i),'Top Right Geodetic Coordinates');
% output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=bound_poly(4,1);
% output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=bound_poly(4,2);

%% Grid
% TODO: Data, as read in provided simulated datasets, does not match predicted frequency bounds
output_meta.Grid.ImagePlane='SLANT';
output_meta.Grid.Type='RGZERO';
output_meta.Grid.Row.Sgn=-1; % Always true for NISAR
output_meta.Grid.Col.Sgn=-1; % Always true for NISAR
output_meta.Grid.Col.KCtr=0;
output_meta.Grid.Row.DeltaKCOAPoly=0;
% Communications with JPL: All SLC data should be uniform weighting
output_meta.Grid.Row.WgtFunct = get_hdf_data(HDF5_fid,...
    '/science/LSAR/SLC/metadata/processingInformation/parameters','rangeChirpWeighting');
if all(output_meta.Grid.Row.WgtFunct(1)==output_meta.Grid.Row.WgtFunct)
    output_meta.Grid.Row.WgtType.WindowName = 'UNIFORM'; % In the sample datasets
else
    output_meta.Grid.Row.WgtType.WindowName = 'UNKNOWN';
end
output_meta.Grid.Col.WgtFunct = get_hdf_data(HDF5_fid,...
    '/science/LSAR/SLC/metadata/processingInformation/parameters','azimuthChirpWeighting');
if all(output_meta.Grid.Col.WgtFunct(1)==output_meta.Grid.Col.WgtFunct)
    output_meta.Grid.Col.WgtType.WindowName = 'UNIFORM'; % In the sample datasets
else
    output_meta.Grid.Col.WgtType.WindowName = 'UNKNOWN';
end
% More subfields added below in "per band" section

%% Timeline
% Using zero Doppler times.  Pulse transmission likely started before this
[collectStart, collectStartFrac]=datenum_w_frac(...
    deblank(get_hdf_data(HDF5_fid, '/science/LSAR/identification', ...
    'zeroDopplerStartTime')'));
[collectEnd, collectEndFrac]=datenum_w_frac(...
    deblank(get_hdf_data(HDF5_fid, '/science/LSAR/identification', ...
    'zeroDopplerEndTime')'));
output_meta.Timeline.CollectStart=collectStart + (collectStartFrac/SECONDS_IN_A_DAY);
output_meta.Timeline.CollectDuration=...
    round((collectEnd-collectStart)*SECONDS_IN_A_DAY) + ... % Convert days to seconds
    (collectEndFrac-collectStartFrac); % Handle fractional seconds
output_meta.Timeline.IPP.Set.TStart=0;
output_meta.Timeline.IPP.Set.TEnd=0; % Apply real value later.  Just a placeholder.
output_meta.Timeline.IPP.Set.IPPStart=uint32(0);
% More subfields added below in "per band" section

%% Position
% Compute polynomial from state vectors
% Get reference time
ref_time = get_ref_time(HDF5_fid,'/science/LSAR/SLC/metadata/orbit/time');
% Times in SICD are with respect to time from start of collect.
ref_time_offset = round((ref_time-collectStart)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
    collectStartFrac; % Handle fractional seconds
state_vector_T  = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/orbit','time')';
state_vector_T = state_vector_T + ref_time_offset; % Make with respect to Timeline.CollectStart
state_vector_pos  = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/orbit','position');
% sv2poly.m shows ways to determine best polynomial order, but we just use 6th here
polyorder=min(6, numel(state_vector_T) - 1);
old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
P_x  = polyfit(state_vector_T, state_vector_pos(1,:), polyorder);
P_y  = polyfit(state_vector_T, state_vector_pos(2,:), polyorder);
P_z  = polyfit(state_vector_T, state_vector_pos(3,:), polyorder);
% We don't use these since they are derivable from the position polynomial
% state_vector_vel = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/orbit','velocity');
% state_vector_acc = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/orbit','acceleration');
% P_vx = polyfit(state_vector_T, state_vector_vel(1,:), polyorder);
% P_vy = polyfit(state_vector_T, state_vector_vel(2,:), polyorder);
% P_vz = polyfit(state_vector_T, state_vector_vel(3,:), polyorder);
% P_ax = polyfit(state_vector_T, state_vector_acc(1,:), polyorder);
% P_ay = polyfit(state_vector_T, state_vector_acc(2,:), polyorder);
% P_az = polyfit(state_vector_T, state_vector_acc(3,:), polyorder);
warning(old_state);
% Store position polynomial
output_meta.Position.ARPPoly.X  = P_x(end:-1:1).';
output_meta.Position.ARPPoly.Y  = P_y(end:-1:1).';
output_meta.Position.ARPPoly.Z  = P_z(end:-1:1).';

%% SCPCOA
% This should almost always be left looking for NISAR.
output_meta.SCPCOA.SideOfTrack=get_hdf_data(HDF5_fid,'/science/LSAR/identification','lookDirection');
output_meta.SCPCOA.SideOfTrack=upper(output_meta.SCPCOA.SideOfTrack(1));
% Most subfields added below in "per band" section.

%% ImageFormation
output_meta.ImageFormation.RcvChanProc=struct('NumChanProc',1,'PRFScaleFactor',1);
output_meta.ImageFormation.ImageFormAlgo='RMA';
output_meta.ImageFormation.TStartProc=0;
output_meta.ImageFormation.TEndProc=output_meta.Timeline.CollectDuration;
output_meta.ImageFormation.STBeamComp='NO';
output_meta.ImageFormation.ImageBeamComp='SV';
output_meta.ImageFormation.AzAutofocus='NO';
output_meta.ImageFormation.RgAutofocus='NO';
% More subfields added below in "per band" section
output_meta.RMA.RMAlgoType='OMEGA_K';
output_meta.RMA.ImageType='INCA';
output_meta.RMA.INCA.DopCentroidCOA=true;
% Get reference time
ref_time = get_ref_time(HDF5_fid,'/science/LSAR/SLC/swaths/zeroDopplerTime');
% Times in SICD are with respect to time from start of collect.
ref_time_offset = round((ref_time-collectStart)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
    collectStartFrac; % Handle fractional seconds
zd_t = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/swaths','zeroDopplerTime') + ref_time_offset;
ss_az_s = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/swaths','zeroDopplerTimeSpacing');
if output_meta.SCPCOA.SideOfTrack(1)=='L'
    zd_t = zd_t(end:-1:1);
    ss_az_s = -ss_az_s;
end
grid_r = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/processingInformation/parameters','slantRange');
% Get reference time
ref_time = get_ref_time(HDF5_fid,'/science/LSAR/SLC/metadata/processingInformation/parameters/zeroDopplerTime');
% Times in SICD are with respect to time from start of collect.
ref_time_offset = round((ref_time-collectStart)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
    collectStartFrac; % Handle fractional seconds
grid_zd = get_hdf_data(HDF5_fid,'/science/LSAR/SLC/metadata/processingInformation/parameters','zeroDopplerTime') + ref_time_offset;

%% Radiometric
cal_fields = {'beta0','gamma0','sigma0'};
sicd_cal_fields = {'BetaZeroSFPoly','GammaZeroSFPoly','SigmaZeroSFPoly'};
cal_data = cell(3,1);  % beta0, gamma0, sigma0
for i = 1:numel(cal_fields)
    cal_data{i} = get_hdf_data(HDF5_fid, '/science/LSAR/SLC/metadata/calibrationInformation/geometry',cal_fields{i}).';
    % Get fill value
    gid = H5G.open(HDF5_fid, '/science/LSAR/SLC/metadata/calibrationInformation/geometry/');
    data = H5D.open(gid, cal_fields{i});
    H5G.close(gid);
    attribute = H5A.open(data, '_FillValue');
    H5D.close(data);
    fill = H5A.read(attribute,H5A.get_type(attribute));
    H5A.close(attribute);
    fill = cast(fill,class(cal_data{i}));
    cal_data{i}(cal_data{i}==fill) = NaN;
    if any(diff(cal_data{i}(isfinite(cal_data{i}))))
        if strcmp(cal_fields{i},'beta0')
            warning('META2SICD_NISAR:unexpectedCalibrationValues',...
                'Beta-0 values expected to be constant.');
        end
    else
        % TODO: Validate conversion factor between NISAR cal fields and
        % SICD cal fields. This assumes Sentinel-1 convention.
        output_meta.Radiometric.(sicd_cal_fields{i}) = 1/mean(cal_data{i}(:),'omitnan')^2;
    end
    % We will derive non-constant Radiometric fields from the constant one
    % (likely beta).
end
output_meta.Radiometric.NoiseLevel.NoiseLevelType = 'ABSOLUTE';

%% Process fields specific to each polarimetric band
band_independent_meta=output_meta; % Values that are consistent across all bands
grouped_meta=cell(sum(cellfun(@(x) size(x,1),pols)),1);
meta_i=0;
for i=1:numel(freqs)
    freq_meta = band_independent_meta;

    freq_id=H5G.open(HDF5_fid,['/science/LSAR/SLC/swaths/frequency' freqs(i)]);
    
    %% Grid
    row_ss_id=H5D.open(freq_id,'slantRangeSpacing');
    freq_meta.Grid.Row.SS=H5D.read(row_ss_id);
    H5D.close(row_ss_id);
    % Col.SS is derived after DRateSFPoly below, rather than used from this
    % given field, so that SICD metadata can be internally consistent:
    % col_ss_id=H5D.open(freq_id,'sceneCenterAlongTrackSpacing');
    % output_meta.Grid.Col.SS=H5D.read(col_ss_id);
    % H5D.close(col_ss_id);
    proc_bw_id=H5D.open(freq_id,'processedRangeBandwidth');
    freq_meta.Grid.Row.ImpRespBW=2*H5D.read(proc_bw_id)/SPEED_OF_LIGHT;
    H5D.close(proc_bw_id);
    freq_meta.Grid.Row.DeltaK1=-freq_meta.Grid.Row.ImpRespBW/2;
    freq_meta.Grid.Row.DeltaK2=-freq_meta.Grid.Row.DeltaK1;

    %% Timeline
    prf_id=H5D.open(freq_id,'nominalAcquisitionPRF');
    prf=H5D.read(prf_id);
    H5D.close(prf_id);
    freq_meta.Timeline.IPP.Set.IPPEnd=uint32(floor(prf*band_independent_meta.Timeline.CollectDuration));
    freq_meta.Timeline.IPP.Set.IPPPoly=[0; prf];
    freq_meta.Timeline.IPP.Set.TEnd=freq_meta.Timeline.CollectDuration;

    %% RadarCollection
    fc_id=H5D.open(freq_id,'acquiredCenterFrequency');
    fc=H5D.read(fc_id);
    H5D.close(fc_id);
    bw_id=H5D.open(freq_id,'acquiredRangeBandwidth');
    bw=H5D.read(bw_id);
    H5D.close(bw_id);
    freq_meta.RadarCollection.TxFrequency.Min=fc-(bw/2);
    freq_meta.RadarCollection.TxFrequency.Max=fc+(bw/2);
    % No waveform parameters provided
    for j=1:size(pols{i},1)
        freq_meta.RadarCollection.RcvChannels.ChanParameters(j).TxRcvPolarization=[pols{i}(j,1) ':' pols{i}(j,2)];
    end
    tx_pol = unique(pols{i}(:,1));
    if numel(tx_pol)==1
        freq_meta.RadarCollection.TxPolarization = tx_pol;
    else
        freq_meta.RadarCollection.TxPolarization = 'SEQUENCE';
        for j = 1:numel(tx_pol)
            freq_meta.RadarCollection.TxSequence.TxStep(j).WFIndex = j;
            freq_meta.RadarCollection.TxSequence.TxStep(j).TxPolarization = tx_pol(j);
        end
    end

    %% ImageFormation
    fc_id=H5D.open(freq_id,'processedCenterFrequency');
    fc=H5D.read(fc_id);
    H5D.close(fc_id);
    bw_id=H5D.open(freq_id,'processedRangeBandwidth');
    bw=H5D.read(bw_id);
    H5D.close(bw_id);
    freq_meta.ImageFormation.TxFrequencyProc.MinProc=fc-(bw/2);
    freq_meta.ImageFormation.TxFrequencyProc.MaxProc=fc+(bw/2);
    range_id=H5D.open(freq_id,'slantRange');
    r_ca_sampled=H5D.read(range_id);
    H5D.close(range_id);
    dop_bw_id=H5D.open(freq_id,'processedAzimuthBandwidth');
    dop_bw=H5D.read(dop_bw_id);
    H5D.close(dop_bw_id);
    H5G.close(freq_id);
    dopcentroid_sampled = get_hdf_data(HDF5_fid,[...
        '/science/LSAR/SLC/metadata/processingInformation/parameters/frequency' freqs(i)], ...
        'dopplerCentroid').';
    doprate_sampled = get_hdf_data(HDF5_fid,[...
        '/science/LSAR/SLC/metadata/processingInformation/parameters/frequency' freqs(i)], ...
        'azimuthFMRate').';

    for j=1:size(pols{i},1)
        pol_meta=freq_meta;

        %% ImageData
        dset_id=H5D.open(HDF5_fid,['/science/LSAR/SLC/swaths/frequency' freqs(i) '/' pols{i}(j,:)]);
        dspace_id=H5D.get_space(dset_id);
        [~, datasize] = H5S.get_simple_extent_dims(dspace_id);
        pol_meta.ImageData.NumCols=uint32(datasize(1));
        pol_meta.ImageData.NumRows=uint32(datasize(2));
        pol_meta.ImageData.FullImage=pol_meta.ImageData;
        if (H5T.get_class(H5D.get_type(dset_id)) == H5ML.get_constant_value('H5T_COMPOUND'))
            if (H5T.get_member_class(H5D.get_type(dset_id),0) == H5ML.get_constant_value('H5T_FLOAT'))
                % Could be either 16- or 32-bit float.  Either way it will
                % cast to this SICD type.
                pol_meta.ImageData.PixelType='RE32F_IM32F';
            elseif (H5T.get_member_class(H5D.get_type(dset_id),0) == H5ML.get_constant_value('H5T_INTEGER')) && ...
                    H5T.get_size(H5D.get_type(dset_id)) == 4
                pol_meta.ImageData.PixelType='RE16I_IM16I';
            else
                error('META2SICD_NISAR:unrecognizedDataType','Must be 16- or 32-bit complex float.');
            end
        else
            error('META2SICD_NISAR:unrecognizedDataType','Must be 16- or 32-bit complex float.');
        end
        H5S.close(dspace_id);
        H5D.close(dset_id);
        pol_meta.ImageData.FirstRow=uint32(0); pol_meta.ImageData.FirstCol=uint32(0);
        % SCP pick is arbitrary.  Just pick something near center.
        pol_meta.ImageData.SCPPixel.Col = floor(datasize(1)/2);
        pol_meta.ImageData.SCPPixel.Row = floor(datasize(2)/2);

        %% ImageFormation
        pol_meta.ImageFormation.RcvChanProc.ChanIndex=meta_i;
        pol_meta.ImageFormation.TxRcvPolarizationProc=...
            freq_meta.RadarCollection.RcvChannels.ChanParameters(j).TxRcvPolarization;
        
        %% RMA
        pol_meta.RMA.INCA.R_CA_SCP = r_ca_sampled(pol_meta.ImageData.SCPPixel.Row+1);
        scp_ca_time = zd_t(pol_meta.ImageData.SCPPixel.Col+1);
        % Compute DRateSFPoly
        % For the purposes of the DRateSFPoly computation, we ignore any
        % changes in velocity or doppler rate over the azimuth dimension.
        pos_coefs = [P_x(:) P_y(:) P_z(:)];
        % Velocity is derivate of position.
        vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
        vel_x = polyval(vel_coefs(:,1), scp_ca_time);
        vel_y = polyval(vel_coefs(:,2), scp_ca_time);
        vel_z = polyval(vel_coefs(:,3), scp_ca_time);
        vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
        r_ca_poly = [pol_meta.RMA.INCA.R_CA_SCP; 1]; % Polynomial representing range as a function of range distance from SCP
        [~, min_ind] = min(abs(grid_zd - scp_ca_time));  % Closest Doppler rate polynomial to SCP
        coords_rg_m = grid_r - pol_meta.RMA.INCA.R_CA_SCP;
        % TODO: Is it the NISAR convention to give the absolute value of
        % Doppler rate? (It should be negative.) Sample datasets appear
        % this way.
        old_warning_state=warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
        dop_rate_poly=polyfit(coords_rg_m, -doprate_sampled(min_ind,:).', 4);
        warning(old_warning_state);
        % figure; plot(r,dop_rate(:,min_ind),r,polyval(dop_rate_poly,r));  % Show closeness of fit
        dop_rate_poly = dop_rate_poly(end:-1:1);
        pol_meta.RMA.INCA.DRateSFPoly = - conv(dop_rate_poly,r_ca_poly) * ... % Multiplication of two polynomials is just a convolution of their coefficients
            SPEED_OF_LIGHT / (2 * fc * vm_ca_sq(1)); % Assumes a SGN of -1
        % Fields dependent on Doppler rate
        pol_meta.Grid.Col.SS = sqrt(vm_ca_sq(1)) * abs(ss_az_s) * pol_meta.RMA.INCA.DRateSFPoly(1,1);
        % Should be close to this:
        % colss_id=H5D.open(freq_id,'sceneCenterAlongTrackSpacing');
        % pol_meta.Grid.Col.SS = H5D.read(colss_id)
        % H5D.close(colss_id);
        pol_meta.Grid.Col.ImpRespBW = ... % Convert to azimuth spatial bandwidth (cycles per meter)
            min(dop_bw*abs(ss_az_s),1)/pol_meta.Grid.Col.SS; % Can't have more bandwidth in data than sample spacing
        pol_meta.RMA.INCA.TimeCAPoly = [scp_ca_time; ... % With respect to start of collect
            ss_az_s/pol_meta.Grid.Col.SS]; % Convert zero doppler spacing from sec/pixels to sec/meters
        % TimeCOAPoly/DopCentroidPoly/DeltaKCOAPoly
        POLY_ORDER = 3; % Order of polynomial which we want to compute
        coords_az_m = (grid_zd - scp_ca_time) * pol_meta.Grid.Col.SS / ss_az_s;
        timeca_sampled = repmat(grid_zd,[1 numel(grid_r)]);
        % TimeCOAPoly=TimeCA+(DopCentroid/dop_rate)
        timecoa_sampled = timeca_sampled+(dopcentroid_sampled./doprate_sampled);
        pol_meta.RMA.INCA.DopCentroidPoly=polyfit2d(dopcentroid_sampled, coords_az_m, coords_rg_m, POLY_ORDER);
        pol_meta.Grid.Col.DeltaKCOAPoly=...
            pol_meta.RMA.INCA.DopCentroidPoly*ss_az_s/pol_meta.Grid.Col.SS;
        pol_meta.Grid.TimeCOAPoly=polyfit2d(timecoa_sampled, coords_az_m, coords_rg_m, POLY_ORDER);

        %% Radiometric
        % In the simulated data (from UAVSAR), we expect these noise values
        % to be only valid for the same incidence angle as NISAR (small
        % part of swath).  We also expect the values to be in dB, not
        % linear.
        nesz = get_hdf_data(HDF5_fid,['/science/LSAR/SLC/metadata/calibrationInformation/frequency' freqs(i) '/' pols{i}(j,:)], 'nes0').';
        % TODO: Validate conversion factor between NISAR cal fields and
        % SICD cal fields. This assumes Sentinel-1 convention (same as
        % Radiometric fields in top part of this code.)
        sigma0sf = 1./(cal_data{3}.^2);
        noise_samples = nesz - (10*log10(sigma0sf));
        pol_meta.Radiometric.NoiseLevel.NoisePoly = ...
            polyfit2d(noise_samples, coords_az_m, coords_rg_m, POLY_ORDER);

        %% GeoData
        % Now that sensor model fields have been populated, we can populate
        % GeoData.SCP more precisely.
        ecf = point_image_to_ground([pol_meta.ImageData.SCPPixel.Row;pol_meta.ImageData.SCPPixel.Col],pol_meta);
        pol_meta.GeoData.SCP.ECF.X=ecf(1);
        pol_meta.GeoData.SCP.ECF.Y=ecf(2);
        pol_meta.GeoData.SCP.ECF.Z=ecf(3);
        llh=ecf_to_geodetic([pol_meta.GeoData.SCP.ECF.X pol_meta.GeoData.SCP.ECF.Y pol_meta.GeoData.SCP.ECF.Z]);
        pol_meta.GeoData.SCP.LLH.Lat=llh(1);
        pol_meta.GeoData.SCP.LLH.Lon=llh(2);
        pol_meta.GeoData.SCP.LLH.HAE=llh(3);

        pol_meta = derived_sicd_fields(pol_meta);
        meta_i=meta_i+1;
        grouped_meta{meta_i}=pol_meta;
    end
end
output_meta = grouped_meta;

end

% Reads HDF5 attribute
function value = get_hdf_attribute(hid_t,path,attribute_name)
gid = H5G.open(hid_t,path);
attribute = H5A.open(gid,attribute_name);
value = H5A.read(attribute,H5A.get_type(attribute));
value = value{1};
H5G.close(gid);
H5A.close(attribute);
end

% Reads HDF5 data
function value = get_hdf_data(hid_t,path,data_name)
gid = H5G.open(hid_t,path);
data = H5D.open(gid,data_name);
value = H5D.read(data);
H5G.close(gid);
H5D.close(data);
end

function ref_time = get_ref_time(hid_t,data_path)
    data = H5D.open(hid_t,data_path);
    attribute = H5A.open(data,'units');
    ref_time_string = H5A.read(attribute,H5A.get_type(attribute))';
    H5A.close(attribute);
    H5D.close(data);
    [ref_time_string,~] = regexp(ref_time_string,'seconds since (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\s*$','tokens','match');
    if numel(ref_time_string)==1  % Should always be the case according to JPL
        ref_time = datenum(ref_time_string{1}{1},'yyyy-mm-dd HH:MM:SS');
    else
        error('META2SICD_NISAR:unrecognizedTimeUnit','Time units must be seconds since a reference UTC time.'); 
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
    if any(isnan(datenum_frac))||isempty(datenum_frac), datenum_frac = 0; end
end

function p = polyfit2d(samples, dim1, dim2, poly_order)
    % Least squares fit for 2D polynomial
    % A*x = b
    [dim1_2d, dim2_2d] = ndgrid(dim1, dim2);
    a = zeros(numel(samples), (poly_order+1)^2);
    for k = 0:poly_order
        for m = 0:poly_order
            a(:,k*(poly_order+1)+m+1) = (dim2_2d(:).^m).*(dim1_2d(:).^k);
        end
    end
    b = zeros((poly_order+1)^2,1);
    for k=1:((poly_order+1)^2)
        b(k)=sum(samples(:).*a(:,k)); % center of aperture
    end
    A=zeros((poly_order+1)^2);
    for k=1:((poly_order+1)^2)
        for m=1:((poly_order+1)^2)
            A(k,m)=sum(a(:,k).*a(:,m));
        end
    end
    old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
    x_coa=A\b; % MATLAB often flags this as badly scaled, but results still appear valid
    warning(old_warning_state);
    p=reshape(x_coa, poly_order+1, poly_order+1);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////