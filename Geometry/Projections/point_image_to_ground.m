function [ gpos ] = point_image_to_ground( im_points, sicd_meta, varargin )
%POINT_IMAGE_TO_GROUND Transforms pixel row, col to ground ECF coordinate
% This function implements the SICD Image Projections Description Document:
% http://www.gwg.nga.mil/ntb/baseline/docs/SICD/index.html
%
% gpos = point_image_to_ground(im_points, sicd_meta, 'PropertyName', PropertyValue, ...)
%
% Inputs:
%    im_points   - [2xN] (row; column) coordinates of N points in image (or
%                  subimage if FirstRow/FirstCol are nonzero).  Zero-based,
%                  following SICD convention (rather than MATLAB
%                  convention, which is one-based); that is, upper-left
%                  pixel is [0;0].
%    sicd_meta   - SICD meta data structure
%
%       Property name     Description
%       projection_type   'plane', 'hae', 'dem'.  Default plane.
%       gref              Ground plane reference point ECF coordinates (m).
%                         Default is SCP.  Only valid if projection_type is
%                         'plane'.
%       ugpn              Ground plane unit normal vector.  Default is
%                         tangent to the surface of constant geodetic
%                         height above the WGS-84 reference ellipsoid
%                         passing through GREF.  Only valid if
%                         projection_type is 'plane'.
%       hae0              Surface height (m) above the WGS-84 reference
%                         ellipsoid for projection point.  Only valid if
%                         projection_type is 'hae'.
%       delta_hae_max     Height threshold for convergence of iterative
%                         constant HAE computation (m).  Default 1 meter.
%                         Only valid if projection_type is 'hae'.
%       hae_nlim          Maximum number of iterations allowed for constant
%                         hae computation.  Default 3.  Only valid if
%                         projection_type is 'hae'.
%       delta_arp         ARP position adjustable parameter (ECF, m).  Default 0.
%       delta_varp        VARP position adjustable parameter (ECF, m/s).  Default 0.
%       range_bias        Range bias adjustable parameter (m).  Default 0.
%       adj_params_frame  Coordinate frame used for expressing delta_arp
%                         and delta_varp adjustable parameters.  Allowed
%                         values: 'ECF', 'RIC_ECF', 'RIC_ECI'. Default ECF.
%
% Outputs:
%    gpos        - [3xN] ECF Ground Points along the R/Rdot contour
%
% Contributors: Thomas McDowall, Harris Corporation
%               Lisa Talbot, NGA/IB
%               Rocco Corsetti, NGA/IB
%               Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
if isvector(im_points), im_points = im_points(:); end; % Assure orientation
scp = [sicd_meta.GeoData.SCP.ECF.X; ...
    sicd_meta.GeoData.SCP.ECF.Y; ...
    sicd_meta.GeoData.SCP.ECF.Z]; % Default values for some input arguments require this
scp_llh = ecf_to_geodetic(scp);
% scp_llh = sicd_meta.GeoData.SCP.LLH; % Doesn't seem to be exactly the same
scp_hae = scp_llh(3);
p = inputParser;
p.addParamValue('projection_type','plane', @(x) any(strcmp(x,{'plane','hae','hae_newton','dem'})));
p.addParamValue('gref',scp, @(x) (size(x,1)==3)||(numel(x)==3)); % ECF (meters)
if isfield(sicd_meta,'PFA')&&isfield(sicd_meta.PFA,'FPN')
    default_ugpn=[sicd_meta.PFA.FPN.X; sicd_meta.PFA.FPN.Y; sicd_meta.PFA.FPN.Z];
else
    default_ugpn=wgs_84_norm(scp);
end
p.addParamValue('ugpn',default_ugpn, @(x) (size(x,1)==3)||(numel(x)==3)); % ECF
p.addParamValue('hae0',scp_hae, @isscalar); % meters
p.addParamValue('delta_hae_max', 1, @(x) isscalar(x) && (x>0)); % meters
p.addParamValue('hae_nlim', 3, @isscalar);
% p.addParamValue('dem',[]); % Not sure how this would work yet, DEM filename?
% Adjustable parameters
p.addParamValue('delta_arp',[0 0 0], @(x) numel(x)==3); % ECF? (meters)
p.addParamValue('delta_varp',[0 0 0], @(x) numel(x)==3); % ECF? (meters/s)
p.addParamValue('range_bias',0, @isscalar); % meters
p.addParamValue('adj_params_frame','ECF', @(x) any(strcmpi(x,{'ECF','RIC_ECF','RIC_ECI'})));
p.FunctionName = mfilename;
p.parse(varargin{:});

%% Chapter 4: Compute projection parameters 
try
    % Computation of r/rdot is specific to image formation type
    [r, rdot, arp_coa, varp_coa] = coa_projection_set(sicd_meta, im_points);
catch % If full metadata not available, try to make the approximation of uniformly sampled grid
    sicd_meta.Grid.Type = 'PLANE';
    if ~isfield(sicd_meta.Grid, 'TimeCOAPoly') % Another approximation that may have to be made
        sicd_meta.Grid.TimeCOAPoly = sicd_meta.Timeline.CollectDuration/2;
    end
    [r, rdot, arp_coa, varp_coa] = coa_projection_set(sicd_meta, im_points);
    warning('point_image_to_ground:IncompleteMetadata',...
        'Unable to compute precise position due to incomplete metadata.  Resorting to approximation.');
end
% After r/rdot is computed, the rest is generic to all SAR.

%% Apply adjustable parameters
if strcmpi(p.Results.adj_params_frame,'ECF')
    % No transformation necessary
    delta_arp = p.Results.delta_arp(:);
    delta_varp = p.Results.delta_varp(:);
else % Translate from RIC frame to ECF frame
    % Use the RIC frame at SCP COA time, not at COA time for im_points
    ARP_SCP_COA = [sicd_meta.SCPCOA.ARPPos.X; ...
        sicd_meta.SCPCOA.ARPPos.Y; ...
        sicd_meta.SCPCOA.ARPPos.Z];
    VARP_SCP_COA = [sicd_meta.SCPCOA.ARPVel.X; ...
        sicd_meta.SCPCOA.ARPVel.Y; ...
        sicd_meta.SCPCOA.ARPVel.Z];
    if strcmpi(p.Results.adj_params_frame,'RIC_ECI')
        T_ECEF_RIC = ric_ecf_mat(ARP_SCP_COA, VARP_SCP_COA, 'eci');
    else % RIC_ECF
        T_ECEF_RIC = ric_ecf_mat(ARP_SCP_COA, VARP_SCP_COA, 'ecf');
    end
    delta_arp = T_ECEF_RIC * p.Results.delta_arp(:);
    delta_varp = T_ECEF_RIC * p.Results.delta_varp(:);
end
arp_coa = arp_coa + repmat(delta_arp,1,size(arp_coa,2));
varp_coa = varp_coa + repmat(delta_varp,1,size(arp_coa,2));
r = r + p.Results.range_bias;

%% Perform actual projection
switch p.Results.projection_type
    case 'plane' % Chapter 5 of SICD document
        gpos = point_to_ground_plane(r, rdot, arp_coa, varp_coa, ...
            p.Results.gref, p.Results.ugpn);
    case 'hae' % Chapter 9 of SICD projection document (draft)
        gpos = point_to_hae(r, rdot, arp_coa, varp_coa, scp, ...
            p.Results.hae0, p.Results.delta_hae_max, p.Results.hae_nlim);
    case 'hae_newton' % Spotlight Synthetic Aperture Radar (SAR) Sensor Model, NGA.SIG.0005_1.0, 2010-03-30
        if ~strcmp(sicd_meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
            error('point_image_to_ground:invalid_projection_mode','HAE_NEWTON projection method only valid for spotlight data.');
        end
        gpos = point_to_hae_newton(r, rdot, arp_coa, varp_coa, scp, ...
            p.Results.hae0);
    case 'dem' % Chapter 10 of SICD document (draft)
        gpos = point_to_DEM(r, rdot, arp_coa, varp_coa, p.Results.dem); % TODO
    otherwise
        % Unrecognized projection type
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////