function [ip, deltaGPn, iter] = point_scene_to_image(S, sicd_meta, varargin)
%POINT_SCENE_TO_IMAGE Transforms a 3D ECF point to pixel row, col
% This function implements the SICD Image Projections Description Document:
% http://www.gwg.nga.mil/ntb/baseline/docs/SICD/index.html
%
% [ip, deltaGPn, iter] = point_scene_to_image(S, sicd_meta, 'PropertyName', PropertyValue, ...)
%
% Inputs:
%    S            - [3xN] N Ground Points in ECF coordinates
%    sicd_meta    - SICD meta data structure
%
%       Property name     Description
%       delta_gp_max      Ground plane displacement tol (m), default = quarter pixel
%       delta_arp         ARP position adjustable parameter (ECF, m).  Default 0.
%       delta_varp        VARP position adjustable parameter (ECF, m/s).  Default 0.
%       range_bias        Range bias adjustable parameter (m).  Default 0.
%       adj_params_frame  Coordinate frame used for expressing delta_arp
%                         and delta_varp adjustable parameters.  Allowed
%                         values: 'ECF', 'RIC_ECF', 'RIC_ECI'. Default ECF.
%
% Outputs:
%    ip           - [2xN] (row; column) coordinates of N points in image
%                   (or subimage if FirstRow/FirstCol are nonzero).
%                   Zero-based, following SICD convention (rather than
%                   MATLAB convention, which is one-based); that is,
%                   upper-left pixel is [0;0].
%    deltaGPn     - Residual ground plane displacement (m)
%    iter         - # iterations required
%
% Authors: Thomas McDowall, Harris Corporation
%          Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Extract the relevant SICD info
% Convert ground / scene point to image coordinates
SCP_Row = double(sicd_meta.ImageData.SCPPixel.Row);
SCP_Col = double(sicd_meta.ImageData.SCPPixel.Col);

First_Row = double(sicd_meta.ImageData.FirstRow);
First_Col = double(sicd_meta.ImageData.FirstCol);
% We override this so that point_image_to_ground will work in coordinate
% space of full image. We will add it back in at the end.
sicd_meta.ImageData.FirstRow = 0;
sicd_meta.ImageData.FirstCol = 0;

Row_SS = sicd_meta.Grid.Row.SS;
Col_SS = sicd_meta.Grid.Col.SS;

% Parse input parameters
p = inputParser;
p.addParamValue('delta_gp_max',0.25 * sqrt(Row_SS^2 + Col_SS^2), @isscalar); % meters
% Adjustable parameters
p.addParamValue('delta_arp',[0 0 0], @(x) numel(x)==3); % ECF? (meters)
p.addParamValue('delta_varp',[0 0 0], @(x) numel(x)==3); % ECF? (meters/s)
p.addParamValue('range_bias',0, @isscalar); % meters
p.addParamValue('adj_params_frame','ECF', @(x) any(strcmpi(x,{'ECF','RIC_ECF','RIC_ECI'})));
p.FunctionName = mfilename;
p.parse(varargin{:});

uRow = [sicd_meta.Grid.Row.UVectECF.X;...
    sicd_meta.Grid.Row.UVectECF.Y;...
    sicd_meta.Grid.Row.UVectECF.Z];

uCol = [sicd_meta.Grid.Col.UVectECF.X;...
    sicd_meta.Grid.Col.UVectECF.Y;...
    sicd_meta.Grid.Col.UVectECF.Z];

%% 3.1 SCP Projection Equations
SCP = [sicd_meta.GeoData.SCP.ECF.X; sicd_meta.GeoData.SCP.ECF.Y; sicd_meta.GeoData.SCP.ECF.Z];

ARP_SCP_COA = [sicd_meta.SCPCOA.ARPPos.X; sicd_meta.SCPCOA.ARPPos.Y; sicd_meta.SCPCOA.ARPPos.Z];
VARP_SCP_COA = [sicd_meta.SCPCOA.ARPVel.X; sicd_meta.SCPCOA.ARPVel.Y; sicd_meta.SCPCOA.ARPVel.Z];

ARP_SCP_COAminusSCP = ARP_SCP_COA - SCP;
%R_SCP_COA = norm(ARP_SCP_COAminusSCP);
%Rdot_SCP_COA = (1/R_SCP_COA) * dot(VARP_SCP_COA, ARP_SCP_COAminusSCP);

% Normal to instantaneous slant plane that contains SCP at SCP COA is
% tangent to R/Rdot contour at SCP. Points away from center of Earth. Use
% LOOK to establish sign.
if strcmp(sicd_meta.SCPCOA.SideOfTrack, 'L')
    LOOK = 1;
else
    LOOK = -1;
end
SPN_SCP_COA = LOOK * cross(VARP_SCP_COA, -ARP_SCP_COAminusSCP);
uSPN = SPN_SCP_COA./repmat(sqrt(sum(SPN_SCP_COA.^2)),3,1);

%% 6.1 Scene To Image: Single Scene Point

% Spherical earth ground plane normal (exact orientation of plane is not
% critical, as long as it contains S)
if isvector(S), S = S(:); end; % Assure orientation
uGPN = S./repmat(sqrt(sum(S.^2)),3,1);

% Ground plane points are projected along straight lines to the image plane. The
% GP to IP direction is along the SCP COA slant plane normal. Also, compute
% image plane unit normal, uIPN. Compute projection scale factor SF.
uPROJ = uSPN;
IPN = cross(uRow, uCol); % should match SICD.PFA.IPN for PFA data
uIPN = IPN./repmat(sqrt(sum(IPN.^2)),3,1); % pointing away from center of earth
SF = dot(uPROJ, uIPN);

% Initialize
Gn = S;

% 2.4 - Image Plane parameters
% The following section is for ground to image with non-orthogonal axes.
% theta col is angle between uRow and uCol if not 0.
cos_theta_col = dot(uRow, uCol);
sin_theta_col = sqrt(1 - cos_theta_col.^2);
GI = (1 / sin_theta_col^2) * [1 -cos_theta_col; -cos_theta_col 1];

done = false;

% Iterate the ground to image transform
iter = ones(1,size(S,2));
to_iter = true(1,size(S,2));
ip = zeros(2,size(S,2));
deltaPn = zeros(3,size(S,2));
deltaGPn = zeros(1,size(S,2));
while(~done && all(iter < 6))
    % Project ground plane point Gn to image plane point In. The projection
    % distance is DISTn. Compute image coordinates xrow and ycol.
    DISTn = (1/SF) * dot((repmat(SCP,1,size(Gn,2)) - Gn), repmat(uIPN,1,size(Gn,2)));
    In = Gn + repmat(DISTn,3,1) .* repmat(uPROJ,1,size(DISTn,2));
    
    % For a point at IPP, corresponding SCP pixel-centered coord is
    deltaIPP = In - repmat(SCP,1,size(In,2));
    ip_iter = GI * ...
        [dot(deltaIPP, repmat(uRow,1,size(In,2))); ...
        dot(deltaIPP, repmat(uCol,1,size(In,2)))];
    
    xrow = ip_iter(1,:);
    ycol = ip_iter(2,:);
    
    irow = xrow / Row_SS;
    icol = ycol / Col_SS;
    
    row = irow + SCP_Row;
    col = icol + SCP_Col;
    
    ip(:,to_iter) = [row; col];
    
    % Transform to ground plane containing the scene point(s) S
    Pn = point_image_to_ground(ip(:,to_iter), sicd_meta, ...
        'projection_type', 'plane', 'gref', S(:,to_iter), 'ugpn', uGPN(:,to_iter), ...
        ... % Pass through adjustable parameters
        'delta_arp', p.Results.delta_arp, ...
        'delta_varp', p.Results.delta_varp, ...
        'range_bias', p.Results.range_bias, ...
        'adj_params_frame', p.Results.adj_params_frame);
    
    %   gndPt = ecfToGeodetic(Pn');
    %   fprintf('Ground Point: %f %f %f\n', gndPt(1), gndPt(2), gndPt(3));
    
    % Compute displacement between Pn and S
    deltaPn(:,to_iter) = S(:,to_iter) - Pn;
    deltaGPn(to_iter) = sqrt(sum(deltaPn(:,to_iter).^2));
    
    old_iter = to_iter;
    to_iter = deltaGPn > p.Results.delta_gp_max; % Need to iterate further on these points
    if any(to_iter)
        Gn = Gn(:,to_iter(old_iter)) + deltaPn(:,to_iter);
        iter(to_iter) = iter(to_iter) + 1;
    else
        done = true;
    end
    
end
if iter > 5
    warning('point_scene_to_image:SolutionDidNotConverge',...
        'Solution for image coordinate did not converge.');
end

ip(1,:) = ip(1,:) - First_Row;
ip(2,:) = ip(2,:) - First_Col;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////