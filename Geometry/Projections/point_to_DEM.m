function [GPP] = point_to_DEM(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, SCP, ...
    delta_HAE_MAX, NLIM, DEM, del_DISTrrc, del_HDlim)
%POINT_TO_DEM transforms center of aperture projection set to DEM in ECF
% coordinates via algorithm in SICD Image Projections.
%
% point_to_DEM(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, SCP, ...
%    delta_HAE_MAX, NLIM, DEMDir, del_DISTrrc, del_HDlim)
%   
% Inputs:
%    R_TGT_COA     - [1xN] range to the ARP at COA
%    Rdot_TGT_COA  - [1xN] range rate relative to the ARP at COA
%    ARP_COA       - [Nx3] aperture reference position at tCOA
%    VARP_COA      - [Nx3] velocity at tCOA
%    SCP           - Scene Center Point
%    delta_HAE_MAX - Height threshold for convergence of iterative
%                    projection sequence.  Recommended value: 1.0 m.
%    NLIM          - Maximum number of HAE iterations allowed.
%                    Recommended value: 3.
%    DEM           - Could be SRTM pathname, or structure with
%                    lats/lons/elevations fields where are all are arrays
%                    of same size and elevation is height above WGS-84
%                    ellipsoid.  Currently no support for UTM or other
%                    coordinate spaces.
%    del_DISTrrc   - Maximum distance between adjacent points along the
%                    R/Rdot contour. Recommended value: 10.0 m.
%    del_HDlim     - Height difference threshold for determining if a point
%                    on the R/Rdot contour is on the DEM surface (m).
%                    Recommended value: 0.001 m.
%
% Outputs:
%    GPP          - [3xN] ECF Ground Plane Point along the R/Rdot contour
%
% Author: Tim Cox, NRL 02082016
%         Wade Schwartzkopf, NGA/R
%         Algorithm Based on "SICD Volume 3, Image Projections Description"
%         NGA.STND.0024-3_2.0 (DRAFT 2014-10-31)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Defaults 
if ~exist('del_DISTrrc','var'), del_DISTrrc = 10; end %meters
if ~exist('del_HDlim','var'), del_HDlim = 0.001; end %meters

% Compute look.  Assume the look for all points is the same.
uGPN = wgs_84_norm(SCP);
LOOK = mode(sign(sum(uGPN.'*cross(ARP_COA - repmat(SCP,1,size(ARP_COA,2)),VARP_COA,1),1)));

if ischar(DEM)  % Assume path to SRTM
    % Load DEM data for 5km buffer around points roughly projected to
    % height above ellipsoid provided in SICD SCP.
    gpos = ecf_to_geodetic(point_to_ground_plane(R_TGT_COA, Rdot_TGT_COA, ...
        ARP_COA, VARP_COA, SCP, uGPN));
    LL = [min(gpos(1,:)) min(gpos(2,:))] - ...
        5*[1 1/cosd(min(gpos(1,:)))]/111.111; % Approximately 5km of buffer
    UR = [max(gpos(1,:)) max(gpos(2,:))] + ...
        5*[1 1/cosd(max(gpos(1,:)))]/111.111; % Approximately 5km of buffer
    DEMstr = DEM; DEM = struct();  % Replace DEM string with structure
    [DEM.lats,DEM.lons,DEM.elevations,DEM.meta] = ...
        get_DEM_heights_region(LL, UR, ...
        'DEMBaseDir', DEMstr, ...
        'correctUndulation', true);  % Return height in HAE rather than MSL
        % Geoid undulation file (with standard naming conversion) should be placed on MATLAB path
    [DEM.lats, DEM.lons] = meshgrid(DEM.lats, DEM.lons);
end

HAEmax = max(DEM.elevations(:));
HAEmin = min(DEM.elevations(:));

% Loop through points
NumPoints = numel(R_TGT_COA);
GPP = zeros(3,NumPoints);

for i=1:NumPoints
    %% Step 1 - Compute the center point and the radius of the R/RDot projection
    %%contour, Rrrc
    VMAG = norm(VARP_COA(:,i));
    uVEL = VARP_COA(:,i)/VMAG;
    cos_DCA = -1*Rdot_TGT_COA(i)/VMAG;
    sin_DCA = sqrt(1-cos_DCA*cos_DCA);
    CTR = ARP_COA(:,i)+R_TGT_COA(i)*cos_DCA*uVEL;
    Rrrc = R_TGT_COA(i)*sin_DCA;

    %% Step 2 - Compute unit vectors uRRX and uRRY
    DEC_ARP = norm(ARP_COA(:,i));
    uUP = ARP_COA(:,i)/DEC_ARP;
    RRY = cross(uUP,uVEL);
    uRRY = RRY./norm(RRY);
    uRRX = cross(uRRY,uVEL);

    %% Step 3 - Project R/Rdot contour to constant height HAEmax
    Aa = point_to_hae(R_TGT_COA(i), Rdot_TGT_COA(i), ARP_COA(:,i), VARP_COA(:,i), SCP, ...
                HAEmax, delta_HAE_MAX, NLIM);
    cos_CAa = dot(Aa-CTR,uRRX)/Rrrc;
    % sin_CAa = LOOK*sqrt(1-cos_CAa*cos_CAa);

    %% Step 4 - Project R/Rdot contour to constant height HAEmin
    Ab = point_to_hae(R_TGT_COA(i), Rdot_TGT_COA(i), ARP_COA(:,i), VARP_COA(:,i), SCP, ...
                HAEmin, delta_HAE_MAX, NLIM);
    cos_CAb = dot(Ab-CTR,uRRX)/Rrrc;
    sin_CAb = LOOK*sqrt(1-cos_CAb*cos_CAb);

    %% Step 5 - Compute the step size for points along R/Rdot contour
    del_cos_RRC = del_DISTrrc*(1/Rrrc)*abs(sin_CAb);
    del_cos_DEM = del_DISTrrc*(1/Rrrc)*(abs(sin_CAb)/cos_CAb);
    del_cos_CA = -1*min([del_cos_RRC del_cos_DEM]);

    %% Step 6 - Compute Number of Points Along R/RDot contour
    NPTS = floor(((cos_CAa - cos_CAb)/del_cos_CA))+2;
        
    %% Step 7 - Compute the set of NPTS along R/RDot contour
    cos_CAn = cos_CAb + (0:(NPTS-1))*del_cos_CA;
    sin_CAn = LOOK*sqrt(1-cos_CAn.*cos_CAn);
    P = CTR + Rrrc*(uRRX*cos_CAn + uRRY*sin_CAn);

    %% Step 8 & 9 - For Each Point convert from ECF to DEM coordinates and compute Delta Height
    llh = ecf_to_geodetic(P);
    h = interp2(DEM.lats,DEM.lons,DEM.elevations,llh(1,:),llh(2,:)); % Is there a smarter interpolation?
    del_h = llh(3,:) - h;

    %% Step 10 - Solve for the points that are on the DEM in increasing height (we may just take one for simplicity)
    % Currently only finds first point (lowest WGS-84 HAE).
    % Finding all solutions would require inspecting all zero crossings.
    close_enough = find(abs(del_h)<del_HDlim,1);
    if ~isempty(close_enough)
        GPP(:,i) = P(:,close_enough);
    else
        zero_cross = find(del_h(1:(end-1)).*del_h(2:end)<1,1);
        if isempty(zero_cross) % No valid solutions
            GPP(:,i) = NaN;
        else
            frac = del_h(zero_cross)/(del_h(zero_cross)-del_h(zero_cross+1));
            cos_CA_S = cos_CAb+(zero_cross-1+frac)*del_cos_CA;    
            sin_CA_S = LOOK*sqrt(1-cos_CA_S*cos_CA_S);
            GPP(:,i) = CTR + Rrrc*(cos_CA_S*uRRX+sin_CA_S*uRRY);
        end
    end
end
end