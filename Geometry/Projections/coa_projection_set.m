function [R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, tCOA] = coa_projection_set(SICD, ip)
% COA_PROJECTION_SET Computes the set of fundamental parameters for
% projecting a pixel down to the ground
%
% [R_TGT_COA, Rdot_TGT_COA, ARP_COA, VARP_COA, tCOA] = coa_projection_set(SICD, ip)
%   
% Inputs:
%    SICD        - SICD meta data structure (see NGA SAR Toolbox read_sicd_meta)
%    ip          - [2xN] (row; column) coordinates of N points in image (or
%                  subimage if FirstRow/FirstCol are nonzero).  Zero-based,
%                  following SICD convention (rather than MATLAB
%                  convention, which is one-based); that is, upper-left
%                  pixel is [0;0].
%
% Outputs:
%    R_TGT_COA    - range to the ARP at COA
%    Rdot_TGT_COA - range rate relative to the ARP at COA
%    ARP_COA      - aperture reference position at tCOA
%    VARP_COA     - velocity at tCOA
%    tCOA         - center of aperture time since CDP start for input ip
%
% Authors: Thomas McDowall, Harris Corporation
%          Rocco Corsetti, NGA/IB
%          Lisa Talbot, NGA/IB
%          Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Check for polar format range az grid
Grid_Type = SICD.Grid.Type;
IFA = SICD.ImageFormation.ImageFormAlgo;

%% Extract the relevant SICD info
SCP_Row = double(SICD.ImageData.SCPPixel.Row);
SCP_Col = double(SICD.ImageData.SCPPixel.Col);

First_Row = double(SICD.ImageData.FirstRow);
First_Col = double(SICD.ImageData.FirstCol);

SCP = [SICD.GeoData.SCP.ECF.X;
       SICD.GeoData.SCP.ECF.Y;
       SICD.GeoData.SCP.ECF.Z];

Row_SS = SICD.Grid.Row.SS;
Col_SS = SICD.Grid.Col.SS;

cT_COA = SICD.Grid.TimeCOAPoly;
M_TCOA = size(cT_COA, 1) - 1; % poly order
N_TCOA = size(cT_COA, 2) - 1; % poly order

cARPx = SICD.Position.ARPPoly.X;
cARPy = SICD.Position.ARPPoly.Y;
cARPz = SICD.Position.ARPPoly.Z;
N_ARP = length(cARPx) - 1; % poly order

%% Convert pixel to meters from SCP
row = double(ip(1,:));
col = double(ip(2,:));

irow = First_Row + row - SCP_Row;
icol = First_Col + col - SCP_Col;

xrow = Row_SS * irow;
ycol = Col_SS * icol;

%% 2.5 Image Grid to COA Parameters

% Compute target pixel time
tCOA = zeros(1,size(ip,2));
for m = 0:M_TCOA
   for n = 0:N_TCOA
      tCOA = tCOA + cT_COA(m+1,n+1) .* (xrow.^m) .* (ycol.^n);
   end
end

% Calculate aperture reference position and velocity at tCOA
ARP_COA  = zeros(3,size(ip,2));
for n = 0:N_ARP
   tCOAn = tCOA.^n;
   ARP_COA(1,:) = ARP_COA(1,:) + cARPx(n+1) * tCOAn;
   ARP_COA(2,:) = ARP_COA(2,:) + cARPy(n+1) * tCOAn;
   ARP_COA(3,:) = ARP_COA(3,:) + cARPz(n+1) * tCOAn;
end

VARP_COA = zeros(3,size(ip,2));
for n = 1:N_ARP
   tCOAnm1 = tCOA.^(n-1);
   VARP_COA(1,:) = VARP_COA(1,:) + n * cARPx(n+1) * tCOAnm1;
   VARP_COA(2,:) = VARP_COA(2,:) + n * cARPy(n+1) * tCOAnm1;
   VARP_COA(3,:) = VARP_COA(3,:) + n * cARPz(n+1) * tCOAnm1;
end

% % Double check for SCP, check Fig 1.1 and 2.7
% if (ip(1) == SCP_Row) && (ip(2) == SCP_Col)
%    
%    disp('SICD.SCPCOA info:')
%    fprintf('tCOA: %f\n', SICD.SCPCOA.SCPTime);
%    fprintf('ARPPos: %f %f %f\n', SICD.SCPCOA.ARPPos.X, SICD.SCPCOA.ARPPos.Y, SICD.SCPCOA.ARPPos.Z);
%    fprintf('ARPVel: %f %f %f\n', SICD.SCPCOA.ARPVel.X, SICD.SCPCOA.ARPVel.Y, SICD.SCPCOA.ARPVel.Z);
%    fprintf('DCA: %f\n', SICD.SCPCOA.DopplerConeAng);
%    fprintf('Slant Range: %f\n', SICD.SCPCOA.SlantRange);
%    fprintf('Radius: %f\n', SICD.SCPCOA.SlantRange * sind(SICD.SCPCOA.DopplerConeAng));
%    fprintf('uRng: %f %f %f\n', SICD.Grid.Row.UVectECF.X, SICD.Grid.Row.UVectECF.Y, SICD.Grid.Row.UVectECF.Z);
%    fprintf('Rdot_COA: %f\n\n', -vec_mag([SICD.SCPCOA.ARPVel.X, SICD.SCPCOA.ARPVel.Y, SICD.SCPCOA.ARPVel.Z]) * cosd(SICD.SCPCOA.DopplerConeAng));
%    
% end
% 

%% 4 Image grid to R/Rdot
if strcmp(Grid_Type, 'RGAZIM')
    % Documentation changes variable names
    rg = xrow;
    az = ycol;

    % Compute range and range rate to the SCP at target pixel COA time
    ARPminusSCP = ARP_COA - repmat(SCP,1,size(ip,2));
    R_SCP_TGT_COA = sqrt(sum(ARPminusSCP.^2));
    Rdot_SCP_TGT_COA = (1./R_SCP_TGT_COA) .* dot(VARP_COA, ARPminusSCP,1);
    
    % fprintf('tCOA: %f\n', tCOA);
    % fprintf('ARPPos: %f %f %f\n', ARP_COA');
    % fprintf('ARPVel: %f %f %f\n', VARP_COA');
    % dca = 180 - angleFrom(unit(VARP_COA), unit(ARPminusSCP)) * 180/pi;
    % fprintf('DCA: %f\n', dca);
    % fprintf('Slant Range: %f\n', R_SCP_TGT_COA);
    % fprintf('Radius: %f\n', R_SCP_TGT_COA * sind(dca));
    % uRng = unit(-ARPminusSCP);
    % fprintf('uRng: %f %f %f\n', uRng);
    % fprintf('Rdot_COA: %f\n', Rdot_SCP_TGT_COA);
    
    %% 4.1 Image grid to R/Rdot for RGAZIM and PFA
    % This computation is dependent on grid type and image formation algorithm.
    if strcmp(IFA, 'PFA')
        % PFA specific SICD metadata
        cPA = SICD.PFA.PolarAngPoly;
        N_PA = length(cPA) - 1; % poly order
        
        cKSF = SICD.PFA.SpatialFreqSFPoly;
        N_KSF = length(cKSF) - 1; % poly order
        
        % Compute polar angle theta and its derivative wrt time at the target pixel COA
        theta_TGT_COA = zeros(1,size(ip,2));
        for n = 0:N_PA
            theta_TGT_COA = theta_TGT_COA + cPA(n+1) * tCOA.^n;
        end
        
        dtheta_dt_TGT_COA = zeros(1,size(ip,2));
        for n = 1:N_PA
            dtheta_dt_TGT_COA = dtheta_dt_TGT_COA + n * cPA(n+1) * tCOA.^(n-1);
        end
        
        % Compute polar aperture scale factor (KSF) and derivative wrt polar angle
        KSF_TGT_COA = zeros(1,size(ip,2));
        for n = 0:N_KSF
            KSF_TGT_COA = KSF_TGT_COA + cKSF(n+1) * theta_TGT_COA.^n;
        end
        
        dKSF_dtheta_TGT_COA = zeros(1,size(ip,2));
        for n = 1:N_KSF
            dKSF_dtheta_TGT_COA = dKSF_dtheta_TGT_COA + n * cKSF(n+1) * theta_TGT_COA.^(n-1);
        end
        
        % Compute spatial frequency domain phase slopes in Ka and Kc directions
        % Note: sign for the phase may be ignored as it is cancelled in a
        % subsequent computation.
        dphi_dKa_TGT_COA =  rg .* cos(theta_TGT_COA) + az .* sin(theta_TGT_COA);
        dphi_dKc_TGT_COA = -rg .* sin(theta_TGT_COA) + az .* cos(theta_TGT_COA);
        
        % Compute range relative to SCP
        deltaR_TGT_COA = KSF_TGT_COA .* dphi_dKa_TGT_COA;
        
        % Compute derivative of range relative to SCP wrt polar angle.
        % Scale by derivative of polar angle wrt time.
        d_deltaR_dtheta_TGT_COA = dKSF_dtheta_TGT_COA .* dphi_dKa_TGT_COA + KSF_TGT_COA .* dphi_dKc_TGT_COA;
        
        delta_Rdot_TGT_COA = d_deltaR_dtheta_TGT_COA .* dtheta_dt_TGT_COA;
    %% 4.2 Image grid to R/Rdot for RGAZIM and RGAZCOMP
    elseif strcmp(IFA, 'RGAZCOMP')
        % RGAZCOMP specific SICD metadata
        AzSF = SICD.RgAzComp.AzSF;
        
        deltaR_TGT_COA = rg;
        delta_Rdot_TGT_COA = - sqrt(sum(VARP_COA.^2)) .* (AzSF * az);
    else
        error('coa_projection_set:UnrecognizedProjection',...
            'Unrecognized grid/IFP type.');
    end
    
    % Compute the range and range rate relative to the ARP at COA. The
    % projection to 3D scene point for grid location (rg, az) is along this
    % R/Rdot contour.
    R_TGT_COA = R_SCP_TGT_COA + deltaR_TGT_COA;
    
    Rdot_TGT_COA = Rdot_SCP_TGT_COA + delta_Rdot_TGT_COA;
%% 4.3 Image grid to R/Rdot for RGZERO
elseif strcmp(Grid_Type, 'RGZERO')
    % RMA specific SICD metadata
    cT_CA = SICD.RMA.INCA.TimeCAPoly;
    N_TCA = size(cT_CA,1)-1;
    
    R_CA_SCP = SICD.RMA.INCA.R_CA_SCP;
    
    cDRSF = SICD.RMA.INCA.DRateSFPoly;
    M_DRSF = size(cDRSF,1)-1;
    N_DRSF = size(cDRSF,2)-1;
    
    % Documentation changes variable names
    rg = xrow;
    az = ycol;
    
    % Compute range at closest approach and time of closest approach
    R_CA_TGT = R_CA_SCP + rg;
    
    t_CA_TGT = zeros(1,size(ip,2));
    for n = 0:N_TCA
        t_CA_TGT = t_CA_TGT + cT_CA(n+1) * (az.^n);
    end
    
    % Compute ARP velocity at t_CA_TGT, and magnitude
    VARP_CA_TGT = zeros(3,size(ip,2));
    for n = 1:N_ARP
        t_CA_TGT_nm1 = t_CA_TGT.^(n-1);
        VARP_CA_TGT(1,:) = VARP_CA_TGT(1,:) + n * cARPx(n+1) * t_CA_TGT_nm1;
        VARP_CA_TGT(2,:) = VARP_CA_TGT(2,:) + n * cARPy(n+1) * t_CA_TGT_nm1;
        VARP_CA_TGT(3,:) = VARP_CA_TGT(3,:) + n * cARPz(n+1) * t_CA_TGT_nm1;
    end
    
    VM_CA_TGT = sqrt(sum(VARP_CA_TGT.^2));
    
    % Compute the Doppler Rate Scale Factor for image grid location
    DRSF_TGT = zeros(1,size(ip,2));
    for  m = 0:M_DRSF
        for n = 0:N_DRSF            
            DRSF_TGT = DRSF_TGT + cDRSF(m+1,n+1) .* (rg.^m) .* (az.^n);            
        end
    end
    
    % Compute time difference between COA time and CA time
    dt_COA_TGT = tCOA - t_CA_TGT;
    
    R_TGT_COA = sqrt(R_CA_TGT.^2 + DRSF_TGT .* (VM_CA_TGT.^2) .* (dt_COA_TGT.^2));
    
    Rdot_TGT_COA = (DRSF_TGT./R_TGT_COA) .* (VM_CA_TGT.^2) .* dt_COA_TGT;
%% 4.4-4.6 Image grid to R/Rdot for uniformly spaced grids
elseif any(strcmp(Grid_Type, {'XRGYCR','XCTYAT','PLANE'}))
    % These grids are all uniformly spaced locations in the image plane.
    % They can all be treated the same.
    % SICD metadata
    uRow = [SICD.Grid.Row.UVectECF.X;...
        SICD.Grid.Row.UVectECF.Y;...
        SICD.Grid.Row.UVectECF.Z];
    uCol = [SICD.Grid.Col.UVectECF.X;...
        SICD.Grid.Col.UVectECF.Y;...
        SICD.Grid.Col.UVectECF.Z];
    
    % Image plane point
    IPP = repmat(SCP,1,size(ip,2)) + (uRow*xrow) + (uCol*ycol);

    % Compute R/Rdot
    ARPminusIPP = ARP_COA - IPP;
    R_TGT_COA = sqrt(sum(ARPminusIPP.^2));
    Rdot_TGT_COA = (1./R_TGT_COA) .* dot(VARP_COA, ARPminusIPP,1);
else
    error('coa_projection_set:UnrecognizedProjection',...
        'Unrecognized grid/image formation type.');
end

% fprintf('\nSlant Range to Target: %f\n', R_TGT_COA);
% fprintf('Rdot Target: %f\n\n', Rdot_TGT_COA);

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////