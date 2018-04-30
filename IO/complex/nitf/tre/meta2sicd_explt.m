function [ output_meta ] = meta2sicd_explt( explt_struct )
%META2SICD_EXPLT Converts metadata stored in the EXPLTA or EXPLTB NITF TREs
% into a SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% We will allow these to be populated by derived_sicd_fields later:
% if ~isnan(explt_struct.GRAZE_ANG)&&explt_struct.GRAZE_ANG~=0
%     output_meta.SCPCOA.GrazeAng=explt_struct.GRAZE_ANG;
%     output_meta.SCPCOA.IncidenceAng=90-explt_struct.GRAZE_ANG;
% end
% if ~isnan(explt_struct.SLOPE_ANG)&&explt_struct.SLOPE_ANG~=0
%     output_meta.SCPCOA.SlopeAng=explt_struct.SLOPE_ANG;
%     output_meta.SCPCOA.TwistAng = -acosd( cosd(explt_struct.SLOPE_ANG) / cosd(explt_struct.GRAZE_ANG) ) * sign(explt_struct.SQUINT_ANGLE);
%     if ~isreal(output_meta.SCPCOA.TwistAng) % cos(slope)>cos(graze)=>acos is imaginary
%         output_meta.SCPCOA=rmfield(output_meta.SCPCOA,'TwistAng');
%     end
% end
output_meta.ImageFormation.TxRcvPolarizationProc=[explt_struct.POLAR(1)...
    ':' explt_struct.POLAR(2)];
switch explt_struct.MODE(2:3)
    case 'SP'
        output_meta.Grid.ImagePlane='SLANT';
    case 'GP'
        output_meta.Grid.ImagePlane='GROUND';
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////