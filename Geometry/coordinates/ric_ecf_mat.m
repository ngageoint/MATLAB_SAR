function [T_ECEF_RIC] = ric_ecf_mat(R_ARP,V_ARP,frame_type)
% RIC_ECF_MAT Compute ECF transformation matrix for RIC frame
%
% Author: Rocco Corsetti, NGA/IB
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if strcmpi(frame_type, 'eci') % RIC_ECI frame
    w_e = 7.292115E-5; % Earth inertial spin rate in radians/second, not including precession
    Vi = V_ARP(:) + cross([0 ; 0 ; w_e], R_ARP(:));
elseif strcmpi(frame_type, 'ecf') % RIC_ECF frame
    Vi = V_ARP(:);
end

Rnum = R_ARP(:);
Rden = norm(Rnum);

Cnum = cross(R_ARP(:), Vi);
Cden = norm(Cnum);

R = Rnum/Rden;
C = Cnum/Cden;
I = cross(C, R);

T_ECEF_RIC = [R I C];

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////