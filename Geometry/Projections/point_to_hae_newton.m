function [SPP] = point_to_hae_newton(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, SCP, HAE0, delta_MAX, NLIM)
% POINT_TO_HAE_NEWTON transforms pixel row, col to a constant height
% above the ellipsoed via algorithm in SICD Image Projections.
%
% GPP = point_to_hae_newton(R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, SCP, HAE0)
%
% Inputs:
%    R_TGT_COA    - range to the ARP at COA
%    Rdot_TGT_COA - range rate relative to the ARP at COA
%    ARP_COA      - aperture reference position at tCOA
%    VARP_COA     - velocity at tCOA
%    SCP          - scene center point (ECF meters)
%    HAE0         - Surface height (m) above the WGS-84 reference ellipsoid
%                   for projection point SPP
%    delta_MAX    - Threshold for convergence of iterative projection
%                   sequence.
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

if ~exist('delta_MAX','var')
    delta_MAX = 1e-3;
end
if ~exist('NLIM','var')
    NLIM = 10;
end
SPP = zeros(3,numel(R_TGT_COA));

%% 9. Precise R/Rdot To Constant HAE Surface Projection

% This way of looping through points really adds no efficiency in MATLAB.
% Really we shouled vectorize the matrix operations and not loop.  However,
% this is an easy way to make its behaviour match the other projection
% functions in the toolbox (point_to_ground_plane, point_to_hae).
for points_i = 1:numel(R_TGT_COA)
    
    % Rename variables from SICD projection document notation to "Spotlight
    % Synthetic Aperture Radar (SAR) Sensor Model" notation
    R = SCP(:);
    h = HAE0;
    
    Range_obs = R_TGT_COA(points_i);
    Doppler_obs = -Rdot_TGT_COA(points_i);
    
    R_ARP_curr = ARP_COA(:,(points_i));
    V_ARP_curr = VARP_COA(:,(points_i));
    
    
    R0 = R;
    a = 6378137.0; % WGS84 semimajor axis (m)
    b = 6356752.314245; % WGS84 semiminor axis (m)
    
    delta = Inf;
    iters = 0;
    while (max(abs(delta)) > delta_MAX) && (iters < NLIM)
        
        F1 = -9999*ones(3,3);
        F2 = -9999*ones(3,3);
        F3 = -9999*ones(3,3);
        
        for i = 1:3
            for j = 1:3
                switch j
                    case(1)
                        eps = -0.5;
                    case(2)
                        eps = 0;
                    case(3)
                        eps = +0.5;
                end
                
                R = R0;
                R(i) = R(i) + eps;
                X = R(1); Y = R(2); Z = R(3);
                
                r_est = R - R_ARP_curr;
                
                Range_est = norm(r_est);
                Doppler_est = dot(V_ARP_curr,r_est)/norm(r_est);
                
                F1(i,j) = Range_obs - Range_est;
                F2(i,j) = Doppler_obs - Doppler_est;
                F3(i,j) = ((X^2 + Y^2)/(a+h)^2) + ((Z^2)/(b+h)^2) - 1;
            end
        end
        
        dF1dX = (F1(1,3) - F1(1,1)) / (2*0.5);
        dF1dY = (F1(2,3) - F1(2,1)) / (2*0.5);
        dF1dZ = (F1(3,3) - F1(3,1)) / (2*0.5);
        
        dF2dX = (F2(1,3) - F2(1,1)) / (2*0.5);
        dF2dY = (F2(2,3) - F2(2,1)) / (2*0.5);
        dF2dZ = (F2(3,3) - F2(3,1)) / (2*0.5);
        
        dF3dX = (F3(1,3) - F3(1,1)) / (2*0.5);
        dF3dY = (F3(2,3) - F3(2,1)) / (2*0.5);
        dF3dZ = (F3(3,3) - F3(3,1)) / (2*0.5);
        
        f = -[F1(1,2)
            F2(1,2)
            F3(1,2)];
        
        B = [dF1dX dF1dY dF1dZ
            dF2dX dF2dY dF2dZ
            dF3dX dF3dY dF3dZ];
        
        delta = B\f;
        
        R0 = R0 + delta;
        
        iters = iters + 1;
    end
    
    SPP(:,points_i) = R0; % Convert back to notation from SICD projection document
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////