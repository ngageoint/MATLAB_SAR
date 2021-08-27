function [sicdmeta] = bp_sicd_meta(meta, nbdata, ifp_params)
%BP_SICD_META Creates a SICD metadata structure for an image created by bp_file.m
%
% Inputs:
%     meta:          Collection narrowband data.  This is the same
%                    structure returned by the get_meta function for a
%                    object returned by open_ph_reader.
%     nbdata:        Per-pulse narrowband data.  This is the same structure
%                    returned as the second return parameter of the
%                    read_cphd method of objects created by open_ph_reader.
%     ifp_params:    Structure of some parameters computed during 
%                    bp_image_formation.m
%
% Outputs:
%     sicdmeta:      Image formation and collection information in a
%                    structure compatible with the SICD writer. 
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Create SICD metadata structure.  This is not really possible, since
% back-projection is not supported by SICD, but we will fill it in where
% possible.
sicdmeta = meta2sicd_cphdx(meta,nbdata,ifp_params.channel);

% Image Creation
sicdmeta.ImageCreation.DateTime = now();
sicdmeta.ImageCreation.Profile  = 'bp_sicd_meta.m';

% ImageData
sicdmeta.ImageData.PixelType               = 'RE32F_IM32F';
if isempty(ifp_params.grid)
    sicdmeta.ImageData.NumRows                 = ifp_params.image_size_pixels(1);
    sicdmeta.ImageData.NumCols                 = ifp_params.image_size_pixels(2);
else
    sicdmeta.ImageData.NumRows                 = size(ifp_params.grid,2);
    sicdmeta.ImageData.NumCols                 = size(ifp_params.grid,1);
end
sicdmeta.ImageData.FirstRow                = 0;
sicdmeta.ImageData.FirstCol                = 0;
sicdmeta.ImageData.FullImage.NumRows       = sicdmeta.ImageData.NumRows;
sicdmeta.ImageData.FullImage.NumCols       = sicdmeta.ImageData.NumCols;
sicdmeta.ImageData.SCPPixel.Row            = floor(sicdmeta.ImageData.NumRows/2);
sicdmeta.ImageData.SCPPixel.Col            = floor(sicdmeta.ImageData.NumCols/2);

% GeoData
sicdmeta.GeoData.EarthModel='WGS_84';
if isempty(ifp_params.grid)
    sicdmeta.GeoData.SCP.ECF.X=mean(ifp_params.center(1));
    sicdmeta.GeoData.SCP.ECF.Y=mean(ifp_params.center(2));
    sicdmeta.GeoData.SCP.ECF.Z=mean(ifp_params.center(3));
    scp_lla = ecf_to_geodetic([sicdmeta.GeoData.SCP.ECF.X sicdmeta.GeoData.SCP.ECF.Y sicdmeta.GeoData.SCP.ECF.Z]);
    sicdmeta.GeoData.SCP.LLH.Lat=scp_lla(1);
    sicdmeta.GeoData.SCP.LLH.Lon=scp_lla(2);
    sicdmeta.GeoData.SCP.LLH.HAE=scp_lla(3);
end
switch ifp_params.grid_type
    case {'ground','slant'}
        corners=zeros(2,2,3);
        [x_vec,y_vec] = ndgrid( ifp_params.image_size_meters(2)*[-1/2, 1/2],...
            ifp_params.image_size_meters(1)*[-1/2, 1/2] );
        for i=1:3
           corners(:,:,i) = ifp_params.center(i) +...
               ifp_params.col_unit_vector(i)*x_vec + ifp_params.row_unit_vector(i)*y_vec;
        end
        FRFC = ecf_to_geodetic(corners(1,1,:));
        sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lat=FRFC(1);
        sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lon=FRFC(2);
        FRLC = ecf_to_geodetic(corners(end,1,:));
        sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lat=FRLC(1);
        sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lon=FRLC(2);
        LRLC = ecf_to_geodetic(corners(end,end,:));
        sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lat=LRLC(1);
        sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lon=LRLC(2);
        LRFC = ecf_to_geodetic(corners(1,end,:));
        sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lat=LRFC(1);
        sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lon=LRFC(2);
    case 'DEM'
        sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lat=ifp_params.row_coords(1);
        sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lon=ifp_params.col_coords(1);
        sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lat=ifp_params.row_coords(1);
        sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lon=ifp_params.col_coords(end);
        sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lat=ifp_params.row_coords(end);
        sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lon=ifp_params.col_coords(end);
        sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lat=ifp_params.row_coords(end);
        sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lon=ifp_params.col_coords(1);
end

% Grid
switch ifp_params.grid_type
    case 'ground'
        sicdmeta.Grid.ImagePlane = 'GROUND';
        sicdmeta.Grid.Type = 'PLANE';
    case 'slant'
        sicdmeta.Grid.ImagePlane = 'SLANT';
        sicdmeta.Grid.Type = 'PLANE';
    otherwise % DEM or arbitrary grid
        % This is what really makes back-projection "un-SICD-able".  SICD
        % doesn't handle formation to anything other than a plane.
        sicdmeta.Grid.ImagePlane = 'OTHER';
        sicdmeta.Grid.Type = 'OTHER';
end
pulse_bandwidth = nbdata.SCSS * (numel(ifp_params.sample_range)-1);
subset_Fx0 = nbdata.SC0 + (nbdata.SCSS * double(ifp_params.sample_range(1)-1));
resolution = pulse_info_to_resolution_extent(...
    nbdata.TxPos - nbdata.SRPPos,... % Line-of-site vectors
    subset_Fx0 + (pulse_bandwidth/2),... % Center frequency
    nbdata.SCSS,... % Frequency step size
    pulse_bandwidth); % Bandwidth
if isempty(ifp_params.grid) % If arbitrary grid is passed, none of the following is valid
    sicdmeta.Grid.Row.SS  = ifp_params.image_size_meters(1)/sicdmeta.ImageData.NumRows;
    sicdmeta.Grid.Col.SS  = ifp_params.image_size_meters(2)/sicdmeta.ImageData.NumCols;
    sicdmeta.Grid.Row.Sgn = -1; % Behaves like this for slant-plane formed images.
    sicdmeta.Grid.Col.Sgn = -1; % Not really sure what this means for back-projected images though.
    sicdmeta.Grid.Row.ImpRespWid = 0.886*resolution(1);
    sicdmeta.Grid.Col.ImpRespWid = 0.886*resolution(2);
    % ImpRespBW will be computed later from ImpRespWid in derived_sicd_fields()
    sicdmeta.Grid.Row.WgtType.WindowName = 'UNIFORM';
    sicdmeta.Grid.Col.WgtType.WindowName = 'UNIFORM';
    if all(isfield(ifp_params,{'row_unit_vector','col_unit_vector'}))
        sicdmeta.Grid.Row.UVectECF=cell2struct(num2cell(ifp_params.row_unit_vector(:)),{'X','Y','Z'});
        sicdmeta.Grid.Col.UVectECF=cell2struct(num2cell(ifp_params.col_unit_vector(:)),{'X','Y','Z'});
    end
end

% ImageFormation
sicdmeta.ImageFormation.RcvChanProc.NumChanProc = 1; % This function only handles single-channel data
sicdmeta.ImageFormation.RcvChanProc.ChanIndex = ifp_params.channel;
sicdmeta.ImageFormation.ImageFormAlgo = 'BACKPROJ'; % Not actually supported in SICD
sicdmeta.ImageFormation.TStartProc = nbdata.TxTime(1);
sicdmeta.ImageFormation.TEndProc = nbdata.TxTime(end);
sicdmeta.ImageFormation.TxFrequencyProc.MinProc = mean(subset_Fx0);
sicdmeta.ImageFormation.TxFrequencyProc.MaxProc = sicdmeta.ImageFormation.TxFrequencyProc.MinProc + pulse_bandwidth;
sicdmeta.ImageFormation.STBeamComp = 'NO';
sicdmeta.ImageFormation.ImageBeamComp = 'NO';
sicdmeta.ImageFormation.AzAutofocus = 'NO';
sicdmeta.ImageFormation.RgAutofocus = 'NO';

% SCPCOA
sicdmeta.SCPCOA.SCPTime = nbdata.TxTime(floor(length(nbdata.TxTime)/2)); % Really only valid for all pixels if spotlight

sicdmeta = derived_sicd_fields(sicdmeta); % Computes many derived fields

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////