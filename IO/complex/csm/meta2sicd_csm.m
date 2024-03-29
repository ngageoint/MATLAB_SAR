function [ output_meta ] = meta2sicd_csm( HDF5_fid )
%META2SICD_CSM Converts CSM HDF5 into a SICD-style metadata structure
%
% Also handles KOMPSAT-5, which is nearly identical format
%
% Takes as input a file identifier to an open HDF5 file-- as is returned by
% H5F.open
%
% SICD fields not currently computed:
%    ValidData
%
% Written by: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

SECONDS_IN_A_DAY = 24*60*60;

mission_id = deblank(get_hdf_attribute(HDF5_fid,'Mission ID')');

%% Setup HDF ids
numbands = H5G.get_num_objs(HDF5_fid);
if strcmp( mission_id, 'CSG')
    numbands = 1;  % CSG only has one SICD per HDF5, even for quadpol
end
for i=1:numbands % "pingpong" mode has multiple polarizations
    groupname=['/S0' num2str(i)];
    group_id(i)=H5G.open(HDF5_fid,groupname);
    switch mission_id
        case {'CSK','KMPS'}
            dset_id(i)=H5D.open(HDF5_fid,[groupname '/SBI']);
        otherwise  % 'CSG'
            dset_id(i)=H5D.open(HDF5_fid,[groupname '/IMG']);
    end
    dspace_id(i)=H5D.get_space(dset_id(i));
end

%% CollectionInfo
output_meta.CollectionInfo.CollectorName=deblank(get_hdf_attribute(HDF5_fid,'Satellite ID')');
[collectStart, collectStartFrac]=datenum_w_frac(...
    get_hdf_attribute(HDF5_fid,'Scene Sensing Start UTC')');
% output_meta.CollectionInfo.CoreName=num2str(get_hdf_attribute(HDF5_fid,'Programmed Image ID'));
output_meta.CollectionInfo.CoreName=[datestr(collectStart, 'ddmmmyy'), ...
    '_', output_meta.CollectionInfo.CollectorName, '_', ...
    datestr(collectStart, 'HHMMSS')];
output_meta.CollectionInfo.CollectType='MONOSTATIC';
acq_mode = deblank(get_hdf_attribute(HDF5_fid,'Acquisition Mode')');
switch mission_id
    case 'CSK'
        switch acq_mode
            case {'HIMAGE','PINGPONG'} % Stripmap
                output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
            case {'WIDEREGION','HUGEREGION'} % ScanSAR
                output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
            case {'ENHANCED SPOTLIGHT','SMART'} % "Spotlight"
                output_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
        end
    case 'KMPS'
        switch acq_mode
            case {'STANDARD','ENHANCED STANDARD'} % Stripmap
                output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
            case {'WIDE SWATH','ENHANCED WIDE SWATH'} % ScanSAR
                output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
            case {'HIGH RESOLUTION', 'ENHANCED HIGH RESOLUTION', 'ULTRA HIGH RESOLUTION'} % "Spotlight"
                output_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
        end
    case 'CSG'
        if strncmp(acq_mode, 'SPOTLIGHT', 9)
            output_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
        elseif any(strcmp(acq_mode, {'STRIPMAP', 'QUADPOL'}))
            output_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
        else
            output_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
        end
end
output_meta.CollectionInfo.RadarMode.ModeID=deblank(get_hdf_attribute(HDF5_fid,'Multi-Beam ID')');
output_meta.CollectionInfo.Classification='UNCLASSIFIED';

%% ImageCreation
output_meta.ImageCreation.DateTime=datenum(get_hdf_attribute(...
    HDF5_fid,'Product Generation UTC')','yyyy-mm-dd HH:MM:SS.FFF');
output_meta.ImageCreation.Profile='Prototype';

%% ImageData
output_meta.ImageData=struct(); % Just a placeholder
% Most subfields added below in "per band" section
% Used for computing SCP later

%% GeoData
if strcmp(deblank(get_hdf_attribute(HDF5_fid,'Ellipsoid Designator')'),'WGS84')
    output_meta.GeoData.EarthModel='WGS_84';
end
% Most subfields added below in "per band" section

%% Grid
if strncmp(get_hdf_attribute(HDF5_fid,'Projection ID')','SLANT',5)
    output_meta.Grid.ImagePlane='SLANT';
    output_meta.Grid.Type='RGZERO';
else
    output_meta.Grid.ImagePlane='GROUND';
end
output_meta.Grid.Row.Sgn=-1; % Always true for CSM
output_meta.Grid.Col.Sgn=-1; % Always true for CSM
fc=get_hdf_attribute(HDF5_fid,'Radar Frequency'); % Center frequency
output_meta.Grid.Row.KCtr=2*fc/SPEED_OF_LIGHT;
output_meta.Grid.Col.KCtr=0;
output_meta.Grid.Row.DeltaKCOAPoly=0;
output_meta.Grid.Row.WgtType.WindowName=upper(deblank(get_hdf_attribute(HDF5_fid,'Range Focusing Weighting Function')'));
if strcmpi(output_meta.Grid.Row.WgtType.WindowName,'HAMMING') % The usual CSM weigting
    rowval = get_hdf_attribute(HDF5_fid,'Range Focusing Weighting Coefficient');
elseif strcmpi(output_meta.Grid.Row.WgtType.WindowName,'GENERAL_COSINE')
    output_meta.Grid.Row.WgtType.WindowName='HAMMING';  % They are effectively the same
    rowval = 1-get_hdf_attribute(HDF5_fid,'Range Focusing Weighting Coefficient');
end
if exist('rowval','var')
    output_meta.Grid.Row.WgtType.Parameter.name = 'COEFFICIENT';
    output_meta.Grid.Row.WgtType.Parameter.value = num2str(rowval);
    output_meta.Grid.Row.WgtFunct = raised_cos_fun(512,rowval);
end
output_meta.Grid.Col.WgtType.WindowName=upper(deblank(get_hdf_attribute(HDF5_fid,'Azimuth Focusing Weighting Function')'));
if strcmpi(output_meta.Grid.Col.WgtType.WindowName,'HAMMING') % The usual CSM weigting
    colval = get_hdf_attribute(HDF5_fid,'Azimuth Focusing Weighting Coefficient');
elseif strcmpi(output_meta.Grid.Col.WgtType.WindowName,'GENERAL_COSINE')
    output_meta.Grid.Col.WgtType.WindowName='HAMMING';  % They are effectively the same
    colval = 1-get_hdf_attribute(HDF5_fid,'Azimuth Focusing Weighting Coefficient');
end
if exist('colval','var')
    output_meta.Grid.Col.WgtType.Parameter.name = 'COEFFICIENT';
    output_meta.Grid.Col.WgtType.Parameter.value = num2str(colval);
    output_meta.Grid.Col.WgtFunct = raised_cos_fun(512,colval);
end
% More subfields added below in "per band" section

%% Timeline
[collectStart, collectStartFrac]=datenum_w_frac(...
    get_hdf_attribute(HDF5_fid,'Scene Sensing Start UTC')');
[collectEnd, collectEndFrac]=datenum_w_frac(...
    get_hdf_attribute(HDF5_fid,'Scene Sensing Stop UTC')');
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
output_meta.Timeline.IPP.Set.TEnd=0; % Apply real value later.  Just a placeholder.
output_meta.Timeline.IPP.Set.IPPStart=uint32(0);
% More subfields added below in "per band" section

%% Position
% Compute polynomial from state vectors
[ref_time, ref_time_frac]=datenum_w_frac(...
    get_hdf_attribute(HDF5_fid,'Reference UTC')');
% Times in SICD are with respect to time from start of collect, but
% time in CSM are generally with respect to reference time.
ref_time_offset = round((ref_time-collectStart)*SECONDS_IN_A_DAY) + ... % Convert from days to secs
    (ref_time_frac-collectStartFrac); % Handle fractional seconds
state_vector_T  = get_hdf_attribute(HDF5_fid,'State Vectors Times').'; % In seconds
state_vector_T = state_vector_T + ref_time_offset; % Make with respect to Timeline.CollectStart
state_vector_pos  = get_hdf_attribute(HDF5_fid,'ECEF Satellite Position');
% sv2poly.m shows ways to determine best polynomial order, but 5th is almost always best
polyorder=min(5, numel(state_vector_T) - 1);
P_x  = polyfit(state_vector_T, state_vector_pos(1,:), polyorder);
P_y  = polyfit(state_vector_T, state_vector_pos(2,:), polyorder);
P_z  = polyfit(state_vector_T, state_vector_pos(3,:), polyorder);
% We don't use these since they are derivable from the position polynomial
% state_vector_vel = get_hdf_attribute(HDF5_fid,'ECEF Satellite Velocity');
% state_vector_acc = get_hdf_attribute(HDF5_fid,'ECEF Satellite Acceleration');
% P_vx = polyfit(state_vector_T, state_vector_vel(1,:), polyorder);
% P_vy = polyfit(state_vector_T, state_vector_vel(2,:), polyorder);
% P_vz = polyfit(state_vector_T, state_vector_vel(3,:), polyorder);
% P_ax = polyfit(state_vector_T, state_vector_acc(1,:), polyorder);
% P_ay = polyfit(state_vector_T, state_vector_acc(2,:), polyorder);
% P_az = polyfit(state_vector_T, state_vector_acc(3,:), polyorder);
% Store position polynomial
output_meta.Position.ARPPoly.X  = P_x(end:-1:1).';
output_meta.Position.ARPPoly.Y  = P_y(end:-1:1).';
output_meta.Position.ARPPoly.Z  = P_z(end:-1:1).';

%% RadarCollection
switch mission_id
    case {'CSK', 'KMPS'}
        pols=cell(numbands,1);
        for i=1:numbands
            pols{i}=deblank(get_hdf_attribute(group_id(i),'Polarisation')');
        end
    case 'CSG'
        if strcmp(acq_mode, 'QUADPOL')
            pols = {'HH', 'HV', 'VH', 'VV'};
        else
            pols = {deblank(get_hdf_attribute(HDF5_fid,'Polarization')')};
        end
end
for i = 1:numel(pols)
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
% Most subfields added below in "per band" section

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
output_meta.RMA.INCA.FreqZero=fc;

%% SCPCOA
output_meta.SCPCOA.SideOfTrack=get_hdf_attribute(HDF5_fid,'Look Side')';
output_meta.SCPCOA.SideOfTrack=upper(output_meta.SCPCOA.SideOfTrack(1));
% Most subfields added below in "per band" section, but we grab this field
% now so we know to flip from CSM's EARLY-LATE column order to SICD's
% view-from-above column order.

%% Process fields specific to each polarimetric band
band_independent_meta=output_meta; % Values that are consistent across all bands
grouped_meta=cell(numbands,1);
for i=1:numbands
    output_meta=band_independent_meta;
    
    %% Doppler polynomials
    % These polynomials are used later to determine RMA.INCA.DopCentroidPoly
    try
        dop_rate_poly_rg=get_hdf_attribute(HDF5_fid,'Doppler Rate vs Range Time Polynomial');
    catch
        dop_rate_poly_rg=get_hdf_attribute(group_id(i),'Doppler Rate vs Range Time Polynomial');
    end
    dop_rate_poly_rg=dop_rate_poly_rg(1:find(dop_rate_poly_rg~=0,1,'last'));  % strip off zero coefficients
    dop_rate_poly_rg=-sign(dop_rate_poly_rg(1))*dop_rate_poly_rg;  % Must be negative (KS5 makes it positive)
    switch mission_id
        case {'CSK', 'KMPS'}
            t_az_ref = get_hdf_attribute(HDF5_fid,'Azimuth Polynomial Reference Time');
            t_rg_ref = get_hdf_attribute(HDF5_fid,'Range Polynomial Reference Time');
            dop_poly_az=get_hdf_attribute(HDF5_fid,'Centroid vs Azimuth Time Polynomial');
            dop_poly_rg=get_hdf_attribute(HDF5_fid,'Centroid vs Range Time Polynomial');
            dop_poly_rg=dop_poly_rg(1:find(dop_poly_rg~=0,1,'last'));  % strip off zero coefficients
        case 'CSG'
            B0001_meta=H5G.open(group_id(i),'B0001');
            az_ref_time_nozd = get_hdf_attribute(group_id(i),'Azimuth Polynomial Reference Time');
            first_time = get_hdf_attribute(B0001_meta,'Azimuth First Time');
            last_time = get_hdf_attribute(B0001_meta,'Azimuth Last Time');
            az_fit_times = linspace(first_time, last_time, 11);

            geom_dop_cent_poly = get_hdf_attribute(group_id(i),'Doppler Centroid vs Azimuth Time Polynomial - RAW');

            try
                dop_rate_poly = get_hdf_attribute(HDF5_fid,'Doppler Rate vs Azimuth Time Polynomial');
            catch
                dop_rate_poly = get_hdf_attribute(group_id(i),'Doppler Rate vs Azimuth Time Polynomial');
            end
            
            centroid_values = polyval(geom_dop_cent_poly(end:-1:1), az_fit_times - az_ref_time_nozd);
            rate_values = polyval(dop_rate_poly(end:-1:1), az_fit_times - az_ref_time_nozd);
            zd_times = az_fit_times - centroid_values ./ rate_values;
            t_az_ref = get_hdf_attribute(group_id(i),'Azimuth Polynomial Reference Time - ZD');
            dop_poly_az = polyfit(zd_times - t_az_ref, centroid_values, 4);
            dop_poly_az = dop_poly_az(end:-1:1)';
            dop_poly_az = dop_poly_az(1:find(dop_poly_az~=0,1,'last')); % strip off zero coefficients

            t_rg_ref = get_hdf_attribute(group_id(i),'Range Polynomial Reference Time');
            dop_poly_rg = get_hdf_attribute(group_id(i),'Doppler Centroid vs Range Time Polynomial');
            dop_poly_rg = dop_poly_rg(1:find(dop_poly_rg~=0,1,'last'));  % strip off zero coefficients
    end

    %% ImageData
    [num_dims datasize] = H5S.get_simple_extent_dims(dspace_id(i)); % All polarizations should be same size
    output_meta.ImageData.NumCols=uint32(datasize(1));
    output_meta.ImageData.NumRows=uint32(datasize(2));
    output_meta.ImageData.FullImage=output_meta.ImageData;
    output_meta.ImageData.FirstRow=uint32(0); output_meta.ImageData.FirstCol=uint32(0);
    switch H5T.get_class(H5D.get_type(dset_id(i)))
        case 0  % integer
            output_meta.ImageData.PixelType='RE16I_IM16I';
        case 1  % float
            output_meta.ImageData.PixelType='RE32F_IM32F';
    end
    % There are many different options for picking the SCP point.  We chose
    % the point that is closest to the reference zero-doppler and range
    % times in the CSM metadata.
    t_az_first=get_hdf_attribute(dset_id(i),'Zero Doppler Azimuth First Time'); % Zero doppler time of first column
    ss_az_s=get_hdf_attribute(dset_id(i),'Line Time Interval'); % Image column spacing in zero doppler time (seconds)
    output_meta.ImageData.SCPPixel.Col=uint32(round((t_az_ref-t_az_first)/ss_az_s) + 1);
    if output_meta.SCPCOA.SideOfTrack=='L'
        % Order of columns in SICD goes in reverse time for left-looking
        ss_az_s=-ss_az_s;
        output_meta.ImageData.SCPPixel.Col=uint32(output_meta.ImageData.NumCols-output_meta.ImageData.SCPPixel.Col-1);
        % First column in SICD is actually last line in CSM terminology
        t_az_first=get_hdf_attribute(dset_id(i),'Zero Doppler Azimuth Last Time');
    end
    t_rg_first=get_hdf_attribute(dset_id(i),'Zero Doppler Range First Time'); % Range time of first row
    if strncmp(get_hdf_attribute(HDF5_fid,'Product Type'),'SCS',3)
        % 'Column Time Interval' does not exist in detected products.
        ss_rg_s=get_hdf_attribute(dset_id(i),'Column Time Interval'); % Row spacing in range time (seconds) 
        output_meta.ImageData.SCPPixel.Row=uint32(round((t_rg_ref-t_rg_first)/ss_rg_s) + 1);
    end
    % How Lockheed seems to pick the SCP:
    output_meta.ImageData.SCPPixel.Col = floor(datasize(1)/2);
    output_meta.ImageData.SCPPixel.Row = ceil(datasize(2)/2)-1;
    
    %% GeoData
    % Initially, we just seed this with a rough value.  Later we will put
    % in something more precise.
    switch mission_id
        case {'CSK', 'KMPS'}
            latlon=get_hdf_attribute(group_id(i),'Centre Geodetic Coordinates');
        case 'CSG'
            latlon=get_hdf_attribute(HDF5_fid,'Scene Centre Geodetic Coordinates');
    end
    output_meta.GeoData.SCP.LLH.Lat=latlon(1);
    output_meta.GeoData.SCP.LLH.Lon=latlon(2);
    % CSM generally gives HAE as zero.  Perhaps we should adjust this to DEM.
    output_meta.GeoData.SCP.LLH.HAE=latlon(3);
    ecf=geodetic_to_ecf([output_meta.GeoData.SCP.LLH.Lat output_meta.GeoData.SCP.LLH.Lon output_meta.GeoData.SCP.LLH.HAE]);
    output_meta.GeoData.SCP.ECF.X=ecf(1);
    output_meta.GeoData.SCP.ECF.Y=ecf(2);
    output_meta.GeoData.SCP.ECF.Z=ecf(3);
    % Calling derived_sicd_fields at the end will populate these fields
    % with the sensor model, so we don't need to do it here.
    % latlon=get_hdf_attribute(dset_id(i),'Top Left Geodetic Coordinates');
    % output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=latlon(1);
    % output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=latlon(2);
    % latlon=get_hdf_attribute(dset_id(i),'Bottom Left Geodetic Coordinates');
    % output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=latlon(1);
    % output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=latlon(2);
    % latlon=get_hdf_attribute(dset_id(i),'Bottom Right Geodetic Coordinates');
    % output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=latlon(1);
    % output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=latlon(2);
    % latlon=get_hdf_attribute(dset_id(i),'Top Right Geodetic Coordinates');
    % output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=latlon(1);
    % output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=latlon(2);
    
    %% Grid
    output_meta.Grid.Row.SS=get_hdf_attribute(dset_id(i),'Column Spacing');
    % output_meta.Grid.Row.SS=get_hdf_attribute(dset_id(i),'Column Time Interval')*SPEED_OF_LIGHT/2; % Exactly equivalent to above
    % Col.SS is derived after DRateSFPoly below, rather than used from this
    % given field, so that SICD metadata can be internally consistent:
    % output_meta.Grid.Col.SS=get_hdf_attribute(dset_id(i),'Line Spacing');
    output_meta.Grid.Row.ImpRespBW=2*get_hdf_attribute(group_id(i),'Range Focusing Bandwidth')/SPEED_OF_LIGHT;
    output_meta.Grid.Row.DeltaK1=-output_meta.Grid.Row.ImpRespBW/2;
    output_meta.Grid.Row.DeltaK2=-output_meta.Grid.Row.DeltaK1;
    if strcmpi(output_meta.Grid.Row.WgtType.WindowName,'HAMMING') % The usual CSM weigting
        % Computation of broadening factor for uniform window is:
        % 2 * fzero(@(x) (sin(pi*x)/(pi*x)) - (1/sqrt(2)), .1)
        a = str2double(output_meta.Grid.Row.WgtType.Parameter.value); % Generalized Hamming window parameter
        row_broadening_factor = 2*fzero(@(x) ...
            a*(sin(pi*x)/(pi*x)) + ((1-a)*(sin(pi*(x-1))/(pi*(x-1)))/2) + ...
            ((1-a)*(sin(pi*(x+1))/(pi*(x+1)))/2) - a/sqrt(2), .1);
        output_meta.Grid.Row.ImpRespWid = row_broadening_factor/output_meta.Grid.Row.ImpRespBW;
    end

    %% Timeline
    prf=get_hdf_attribute(group_id(i),'PRF');
    output_meta.Timeline.IPP.Set.IPPEnd=uint32(floor(prf*band_independent_meta.Timeline.CollectDuration));
    output_meta.Timeline.IPP.Set.IPPPoly=[0; prf];
    output_meta.Timeline.IPP.Set.TEnd=output_meta.Timeline.CollectDuration;
    
    %% RadarCollection
    % output_meta.RadarCollection.RefFreqIndex=uint32(0); % Absence of this field means all frequencies are true values
    chirp_length=get_hdf_attribute(group_id(i),'Range Chirp Length');
    chirp_rate=abs(get_hdf_attribute(group_id(i),'Range Chirp Rate'));
    bw=chirp_length*chirp_rate;
    output_meta.RadarCollection.TxFrequency.Min=fc-(bw/2);
    output_meta.RadarCollection.TxFrequency.Max=fc+(bw/2);
    output_meta.RadarCollection.Waveform.WFParameters.TxPulseLength=chirp_length;
    output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth=bw;
    output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart=...
        output_meta.RadarCollection.TxFrequency.Min;
    output_meta.RadarCollection.Waveform.WFParameters.TxFMRate=chirp_rate;
    sample_rate=get_hdf_attribute(group_id(i),'Sampling Rate');
    try
        if isnan(get_hdf_attribute(group_id(i),'Reference Dechirping Time'))
            output_meta.RadarCollection.Waveform.WFParameters.RcvDemodType='CHIRP';
            output_meta.RadarCollection.Waveform.WFParameters.RcvFMRate=0;
        else
            output_meta.RadarCollection.Waveform.WFParameters.RcvDemodType='STRETCH';
            output_meta.RadarCollection.Waveform.WFParameters.RcvFMRate=chirp_rate;
        end
    catch
        output_meta.RadarCollection.Waveform.WFParameters.RcvDemodType='CHIRP';
        output_meta.RadarCollection.Waveform.WFParameters.RcvFMRate=0;
    end
    output_meta.RadarCollection.Waveform.WFParameters.RcvWindowLength=...
        double(get_hdf_attribute(group_id(i),'Echo Sampling Window Length'))/sample_rate;
    output_meta.RadarCollection.Waveform.WFParameters.ADCSampleRate=sample_rate;
    
    %% ImageFormation
    output_meta.ImageFormation.RcvChanProc.ChanIndex=i;
    output_meta.ImageFormation.TxFrequencyProc.MinProc=...
        output_meta.RadarCollection.TxFrequency.Min;
    output_meta.ImageFormation.TxFrequencyProc.MaxProc=...
        output_meta.RadarCollection.TxFrequency.Max;
    if strcmp(mission_id, 'CSG')
        band_pol = deblank(get_hdf_attribute(HDF5_fid,'Polarization')');
        output_meta.ImageFormation.TxRcvPolarizationProc = [band_pol(1) ':' band_pol(2)];
    else
        output_meta.ImageFormation.TxRcvPolarizationProc=...
            band_independent_meta.RadarCollection.RcvChannels.ChanParameters(i).TxRcvPolarization;
    end
    
    %% RMA
    t_rg_scp = t_rg_first + (ss_rg_s*double(output_meta.ImageData.SCPPixel.Row)); % Range time to SCP
    output_meta.RMA.INCA.R_CA_SCP = t_rg_scp*SPEED_OF_LIGHT/2;
    t_az_scp = t_az_first + ...
        (ss_az_s*double(output_meta.ImageData.SCPPixel.Col)); % Zero doppler time of SCP
    % Compute DRateSFPoly
    % For the purposes of the DRateSFPoly computation, we ignore any
    % changes in velocity or doppler rate over the azimuth dimension.
    pos_coefs = [P_x(:) P_y(:) P_z(:)];
    % Velocity is derivate of position.
    vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
    scp_ca_time = t_az_scp + ref_time_offset;
    vel_x = polyval(vel_coefs(:,1), scp_ca_time);
    vel_y = polyval(vel_coefs(:,2), scp_ca_time);
    vel_z = polyval(vel_coefs(:,3), scp_ca_time);
    vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
    r_ca = [output_meta.RMA.INCA.R_CA_SCP; 1]; % Polynomial representing range as a function of range distance from SCP
    if exist('ss_rg_s','var')
        dop_rate_poly_rg_shifted=polyshift(dop_rate_poly_rg, t_rg_scp-t_rg_ref);
        dop_rate_poly_rg_scaled=dop_rate_poly_rg_shifted.*...
            (ss_rg_s/output_meta.Grid.Row.SS).^(0:(length(dop_rate_poly_rg)-1)).';
    end
    output_meta.RMA.INCA.DRateSFPoly = - conv(dop_rate_poly_rg_scaled,r_ca) * ... % Multiplication of two polynomials is just a convolution of their coefficients
        SPEED_OF_LIGHT / (2 * fc * vm_ca_sq(1)); % Assumes a SGN of -1
    % Fields dependent on Doppler rate
    output_meta.Grid.Col.SS = sqrt(vm_ca_sq(1)) * abs(ss_az_s) * output_meta.RMA.INCA.DRateSFPoly(1,1);
    % 'Azimuth Focusing Bandwidth' sometimes represents something larger
    % than ImpRespBW for S2 CSK mode.  Not sure what "transition bandwidth"
    % means, with respect to SAR, but it appears to correlate with non-zero
    % Doppler frequency content in SLC.
    switch mission_id
        case 'CSK'
            dop_bw = get_hdf_attribute(group_id(i),'Azimuth Focusing Transition Bandwidth'); % Doppler frequency
        otherwise
            dop_bw = get_hdf_attribute(group_id(i),'Azimuth Focusing Bandwidth'); % Doppler frequency
    end
    output_meta.Grid.Col.ImpRespBW = ... % Convert to azimuth spatial bandwidth (cycles per meter)
        min(dop_bw*abs(ss_az_s),1)/output_meta.Grid.Col.SS; % Can't have more bandwidth in data than sample spacing
    if strcmpi(output_meta.Grid.Col.WgtType.WindowName,'HAMMING') % The usual CSM weigting
        a = str2double(output_meta.Grid.Col.WgtType.Parameter.value); % Generalized Hamming window parameter
        col_broadening_factor = 2*fzero(@(x) ...
            a*(sin(pi*x)/(pi*x)) + ((1-a)*(sin(pi*(x-1))/(pi*(x-1)))/2) + ...
            ((1-a)*(sin(pi*(x+1))/(pi*(x+1)))/2) - a/sqrt(2), .1);
        output_meta.Grid.Col.ImpRespWid = col_broadening_factor/output_meta.Grid.Col.ImpRespBW;
    end
    output_meta.RMA.INCA.TimeCAPoly = [scp_ca_time; ... % With respect to start of collect
        ss_az_s/output_meta.Grid.Col.SS]; % Convert zero doppler spacing from sec/pixels to sec/meters
    % Compute DopCentroidPoly/DeltaKCOAPoly
    if exist('ss_rg_s','var')
        output_meta.RMA.INCA.DopCentroidPoly=zeros(length(dop_poly_rg),length(dop_poly_az));
        % Compute doppler centroid value at SCP
        output_meta.RMA.INCA.DopCentroidPoly(1) = ...
            polyval(dop_poly_rg(end:-1:1), t_rg_scp-t_rg_ref) + ...
            polyval(dop_poly_az(end:-1:1), t_az_scp-t_az_ref) - ...
            mean([dop_poly_az(1) dop_poly_rg(1)]); % These should be identical
        % Shift 1D polynomials to account for SCP
        dop_poly_az_shifted=polyshift(dop_poly_az, t_az_scp-t_az_ref);
        dop_poly_rg_shifted=polyshift(dop_poly_rg, t_rg_scp-t_rg_ref);
        % Scale 1D polynomials to from Hz/s^n to Hz/m^n
        dop_poly_az_scaled=dop_poly_az_shifted.'.*...
            (ss_az_s/output_meta.Grid.Col.SS).^(0:(length(dop_poly_az)-1));
        dop_poly_rg_scaled=dop_poly_rg_shifted.*...
            (ss_rg_s/output_meta.Grid.Row.SS).^(0:(length(dop_poly_rg)-1)).';
        output_meta.RMA.INCA.DopCentroidPoly(2:end,1)=dop_poly_rg_scaled(2:end);
        output_meta.RMA.INCA.DopCentroidPoly(1,2:end)=dop_poly_az_scaled(2:end);
        output_meta.RMA.INCA.DopCentroidCOA=true;
        output_meta.Grid.Col.DeltaKCOAPoly=...
            output_meta.RMA.INCA.DopCentroidPoly*ss_az_s/output_meta.Grid.Col.SS;
    end
    % TimeCOAPoly
    % TimeCOAPoly=TimeCA+(DopCentroid/dop_rate)
    % Since we can't evaluate this equation analytically, we will evaluate
    % samples of it across our image and fit a 2D polynomial to it.
    POLY_ORDER = 2; % Order of polynomial which we want to compute
    grid_samples = POLY_ORDER + 1; % in each dimension
    coords_az_m = linspace(-double(output_meta.ImageData.SCPPixel.Col),...
        double(output_meta.ImageData.NumCols-output_meta.ImageData.SCPPixel.Col-1), grid_samples) * ...
        output_meta.Grid.Col.SS;
    coords_rg_m = linspace(-double(output_meta.ImageData.SCPPixel.Row),...
        double(output_meta.ImageData.NumRows-output_meta.ImageData.SCPPixel.Row-1), grid_samples) * ...
        output_meta.Grid.Row.SS;
    timeca_sampled = sicd_polyval2d(output_meta.RMA.INCA.TimeCAPoly(:).',coords_az_m,coords_rg_m);
    dopcentroid_sampled = sicd_polyval2d(output_meta.RMA.INCA.DopCentroidPoly,coords_az_m,coords_rg_m);
    doprate_sampled = sicd_polyval2d(dop_rate_poly_rg_scaled.',coords_az_m,coords_rg_m);
    timecoapoly_sampled = timeca_sampled+(dopcentroid_sampled./doprate_sampled);
    % Least squares fit for 2D polynomial
    % A*x = b
    [coords_az_m, coords_rg_m] = ndgrid(coords_az_m, coords_rg_m);
    a = zeros(grid_samples^2, (POLY_ORDER+1)^2);
    for k = 0:POLY_ORDER
        for j = 0:POLY_ORDER
            a(:,k*(POLY_ORDER+1)+j+1) = (coords_rg_m(:).^j).*(coords_az_m(:).^k);
        end
    end
    b_coa = zeros((POLY_ORDER+1)^2,1);
    for k=1:((POLY_ORDER+1)^2)
       b_coa(k)=sum(timecoapoly_sampled(:).*a(:,k)); % center of aperture
    end
    A=zeros((POLY_ORDER+1)^2);
    for k=1:((POLY_ORDER+1)^2)
        for j=1:((POLY_ORDER+1)^2)
            A(k,j)=sum(a(:,k).*a(:,j));
        end
    end
    old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
    x=A\b_coa; % MATLAB often flags this as badly scaled, but results still appear valid
    warning(old_warning_state);
    output_meta.Grid.TimeCOAPoly=reshape(x, POLY_ORDER+1, POLY_ORDER+1);

    %% Radiometric
    if strcmp(mission_id,'CSK')  % Unsure how to handle for CSG/KMPS
        if ~strncmpi(get_hdf_attribute(HDF5_fid,'Range Spreading Loss Compensation Geometry')', 'NONE', 4)
            fact = get_hdf_attribute(HDF5_fid,'Reference Slant Range') ^ ...
                (2 * get_hdf_attribute(HDF5_fid,'Reference Slant Range Exponent'));
            % This code commented out below converts from beta_0 to sigma_0,
            % but we will just populate the more natural beta_0 directly, and
            % let derived_sicd_fields compute the rest of the rcs, sigma_0,
            % gamma_0 terms in a sensor independent way later.
            % if ~strncmpi(get_hdf_attribute(HDF5_fid,'Incidence Angle Compensation Geometry')', 'NONE', 4)
            %     fact = fact * sind(get_hdf_attribute(HDF5_fid,'Reference Incidence Angle'));
            % end
            if isequal(get_hdf_attribute(HDF5_fid,'Calibration Constant Compensation Flag'), 0)
                fact = fact * (1 / (get_hdf_attribute(HDF5_fid,'Rescaling Factor')^2));
                fact = fact / get_hdf_attribute(group_id(i),'Calibration Constant');
                output_meta.Radiometric.BetaZeroSFPoly = fact;
            end
        end
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
    
    grouped_meta{i}=output_meta;
end
output_meta=grouped_meta;

%% Close HDF ids
for i=1:numbands
    H5S.close(dspace_id(i));
    H5D.close(dset_id(i));
    H5G.close(group_id(i));
end

end

% Reads and HDF5 attribute
function value = get_hdf_attribute(hid_t,attribute_name)
attribute=H5A.open_name(hid_t,attribute_name);
value = H5A.read(attribute,H5A.get_type(attribute));
H5A.close(attribute);
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
    datenum_s = datenum(datestring,'yyyy-mm-dd HH:MM:SS');
    datenum_frac = str2double(regexp(datestring,'\.\d*','match'));
    if isnan(datenum_frac), datenum_frac = 0; end;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////