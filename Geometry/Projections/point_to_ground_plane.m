function [GPP] = point_to_ground_plane(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, GREF, GPN)
% POINT_TO_GROUND_PLANE transforms pixel row, col to ground plane ECF coordinate
% via algorithm in SICD Image Projections.
%
% GPP = point_to_ground_plane(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, GREF, GPN)
%   
% Inputs:
%    R_TGT_COA    - [1xN] range to the ARP at COA
%    Rdot_TGT_COA - [1xN] range rate relative to the ARP at COA
%    ARP_COA      - [Nx3] aperture reference position at tCOA
%    VARP_COA     - [Nx3] velocity at tCOA
%    GREF         - [3x1] reference point in the plane to which we are projecting
%    GPN          - [3x1] vector normal to the plane to which we are projecting
%
% Outputs:
%    GPP          - [3xN] ECF Ground Plane Point along the R/Rdot contour
%
% Authors: Thomas McDowall, Harris Corporation
%          Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Ground plane description could be of single plane for all points
% (probably the more typical case) or a plane for each point.
if numel(GPN)==3
    GPN = repmat(GPN,1,numel(R_TGT_COA));
end
if numel(GREF)==3
    GREF = repmat(GREF,1,numel(R_TGT_COA));
end

%% 5. Precise R/Rdot to Ground Plane Projection
% Solve for the intersection of a R/Rdot contour and a ground plane.

% unit normal
uZ = GPN./repmat(sqrt(sum(GPN.^2)),3,1);

% ARP distance from plane
ARPz = dot(ARP_COA - GREF, uZ, 1);
ARPz(ARPz>R_TGT_COA) = NaN; % No solution

% ARP ground plane nadir
AGPN = ARP_COA - repmat(ARPz,3,1) .* uZ;

% Compute ground plane distance (G) from ARP nadir to circle of const range
G = sqrt(R_TGT_COA.^2 - ARPz.^2);

% Compute sine and cosine of grazing angle
cos_GRAZ = G ./ R_TGT_COA;
sin_GRAZ = ARPz ./ R_TGT_COA;

%fprintf('Target GrazeAng: %f\n', acos(cos_GRAZ) * 180/pi);

% Velocity components normal to ground plane and parallel to ground plane.
VMag = sqrt(sum(VARP_COA.^2));
Vz = dot(VARP_COA, uZ, 1);
Vx = sqrt(VMag.^2 - Vz.^2); % Note: For Vx = 0, no Solution

% Orient X such that Vx > 0 and compute unit vectors uX and uY
uX = (1./repmat(Vx,3,1)) .* (VARP_COA - repmat(Vz,3,1) .* uZ);
uY = cross(uZ, uX, 1);

% Compute cosine of azimuth angle to ground plane point
cos_AZ = (-Rdot_TGT_COA + Vz .* sin_GRAZ) ./ (Vx .* cos_GRAZ);
cos_AZ(cos_AZ>1|cos_AZ<-1) = NaN; % R/Rdot combination not possible in given plane
% Don't use to get AZ. Sign may be adjusted in next step.
%fprintf('Illumination Azimuth: %f\n', acos(cos_AZ) * 180/pi);

% Compute sine of azimuth angle. Use LOOK to establish sign.
LOOK = sign(sum(GPN.*cross(ARP_COA - GREF,VARP_COA,1)));
sin_AZ = LOOK .* sqrt(1 - cos_AZ.^2);

%fprintf('Illumination Azimuth: %f\n\n', asin(sin_AZ) * 180/pi);

% Compute Ground Plane Point in ground plane and along the R/Rdot contour
GPP = AGPN + (repmat(G .* cos_AZ,3,1) .* uX) + (repmat(G .* sin_AZ,3,1) .* uY); % in ECF

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////