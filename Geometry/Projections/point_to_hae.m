function [SPP] = point_to_hae(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, SCP, HAE0, delta_HAE_MAX, NLIM)
% POINT_TO_HAE transforms pixel row, col to a constant height
% above the ellipsoid via algorithm in SICD Image Projections.
%
% SPP = point_to_hae(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, SCP, HAE0, delta_HAE_MAX, NLIM)
%   
% Inputs:
%    R_TGT_COA    - range to the ARP at COA
%    Rdot_TGT_COA - range rate relative to the ARP at COA
%    ARP_COA      - aperture reference position at tCOA
%    VARP_COA     - velocity at tCOA
%    SCP          - scene center point (ECF meters)
%    HAE0         - Surface height (m) above the WGS-84 reference ellipsoid
%                   for projection point SPP
%    delta_HAE_MAX- Height threshold for convergence of iterative
%                   projection sequence.
%    NLIM         - Maximum number of iterations allowed.
%
% Outputs:
%    SPP          - [3xN] Surface Projection Point position on the HAE0
%                   surface and along the R/Rdot contour
%
% Authors: Rocco Corsetti, NGA/IB
%          Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% 9. Precise R/Rdot To Constant HAE Surface Projection

numpts = numel(R_TGT_COA);

iters = 0;
delta_HAE = Inf;

% (1) Compute the geodetic ground plane normal at the SCP.
uGPN = wgs_84_norm(SCP);
LOOK = sign(sum(uGPN.'*cross(ARP_COA - repmat(SCP,1,size(ARP_COA,2)),VARP_COA,1),1));
GPP_LLH = ecf_to_geodetic(SCP);
% GREF = repmat(SCP,1,numpts) - repmat((repmat(GPP_LLH(3,:),1,numpts) - HAE0),3,1) .* repmat(uGPN,1,numpts);
GREF = bsxfun(@minus,SCP,bsxfun(@times,bsxfun(@minus, GPP_LLH(3,:), HAE0), uGPN));
while all(abs(delta_HAE) > delta_HAE_MAX) && (iters <= NLIM)
    % (2) Compute the precise projection along the R/Rdot contour to Ground Plane n.
    GPP = point_to_ground_plane(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, GREF, uGPN);
    % (3) Compute the unit vector in the increasing height direction
    uGPN = wgs_84_norm(GPP);
    GPP_LLH = ecf_to_geodetic(GPP);
    delta_HAE = GPP_LLH(3,:) - HAE0;
    GREF = GPP - repmat(delta_HAE,3,1) .* uGPN;
    iters = iters + 1;
    % (4) Test for delta_HAE_MAX and NLIM
end
% (5) Compute the unit slant plane normal vector, uSPN, that is tangent to
% the R/Rdot contour at point GPP
SPN = repmat(LOOK,3,1).* cross(VARP_COA, (GPP - ARP_COA));
uSPN = SPN./repmat(sqrt(sum(SPN.^2)),3,1);
% (6) For the final straight line projection, project from point GPP along
% the slant plane normal (as opposed to the ground plane normal that was
% used in the iteration) to point SLP.
SF = sum(uGPN.*uSPN);
SLP = GPP - repmat(delta_HAE ./ SF,3,1) .* uSPN;
% (7) Assign surface point SPP position by adjusting the HAE to be on the
% HAE0 surface.
SPP_LLH = ecf_to_geodetic(SLP);
SPP_LLH(3,:) = HAE0;
SPP = geodetic_to_ecf(SPP_LLH);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////