function [ output_meta ] = meta2sicd_mpdsra( mpdsra_struct, ORP_ECF )
%META2SICD_MPDSRA Converts metadata stored in the MPDSRA NITF TREs
% into a SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Sometimes a valid ORP is not given in this TRE, but is in another TRE in
% the NITF.  Allow for this possibility.
if isfinite(mpdsra_struct.ORP_X)&&mpdsra_struct.ORP_X~=0&&...
        isfinite(mpdsra_struct.ORP_Y)&&mpdsra_struct.ORP_Y&&...
        isfinite(mpdsra_struct.ORP_Z)&&mpdsra_struct.ORP_Z
    output_meta.GeoData.SCP.ECF.X=mpdsra_struct.ORP_X*FEET_TO_METERS;
    output_meta.GeoData.SCP.ECF.Y=mpdsra_struct.ORP_Y*FEET_TO_METERS;
    output_meta.GeoData.SCP.ECF.Z=mpdsra_struct.ORP_Z*FEET_TO_METERS;
    llh=ecf_to_geodetic([mpdsra_struct.ORP_X mpdsra_struct.ORP_Y mpdsra_struct.ORP_Z]*FEET_TO_METERS);
    output_meta.GeoData.SCP.LLH.Lat=llh(1);
    output_meta.GeoData.SCP.LLH.Lon=llh(2);
    output_meta.GeoData.SCP.LLH.HAE=llh(3);
    output_meta.Position.GRPPoly = output_meta.GeoData.SCP.ECF;
    ORP_ECF=[mpdsra_struct.ORP_X mpdsra_struct.ORP_Y mpdsra_struct.ORP_Z]*FEET_TO_METERS;
end
   
% Use valid ORP from this TRE if given.  Otherwise, use ORP passed as
% function parameter.
if exist('ORP_ECF','var')
    ARP_POS_NED=([mpdsra_struct.ARP_POS_N mpdsra_struct.ARP_POS_E mpdsra_struct.ARP_POS_D].')*FEET_TO_METERS;
    if all(~isnan(ARP_POS_NED))&&all(ARP_POS_NED~=0)
        ARP_POS=ned_to_ecf(ARP_POS_NED,ORP_ECF,true);
        output_meta.SCPCOA.ARPPos.X=ARP_POS(1);
        output_meta.SCPCOA.ARPPos.Y=ARP_POS(2);
        output_meta.SCPCOA.ARPPos.Z=ARP_POS(3);
    end
    ARP_VEL_NED=([mpdsra_struct.ARP_VEL_N mpdsra_struct.ARP_VEL_E mpdsra_struct.ARP_VEL_D].')*FEET_TO_METERS;
    if all(~isnan(ARP_VEL_NED))&&all(ARP_VEL_NED~=0)
        ARP_VEL=ned_to_ecf(ARP_VEL_NED,ORP_ECF,false);
        output_meta.SCPCOA.ARPVel.X=ARP_VEL(1);
        output_meta.SCPCOA.ARPVel.Y=ARP_VEL(2);
        output_meta.SCPCOA.ARPVel.Z=ARP_VEL(3);
    end
    ARP_ACC_NED=([mpdsra_struct.ARP_ACC_N mpdsra_struct.ARP_ACC_E mpdsra_struct.ARP_ACC_D].')*FEET_TO_METERS;
    if all(~isnan(ARP_ACC_NED))&&all(ARP_ACC_NED~=0)
        ARP_ACC=ned_to_ecf(ARP_ACC_NED,ORP_ECF,false);
        output_meta.SCPCOA.ARPAcc.X=ARP_ACC(1);
        output_meta.SCPCOA.ARPAcc.Y=ARP_ACC(2);
        output_meta.SCPCOA.ARPAcc.Z=ARP_ACC(3);
    end
end
% MPDSRA ORP pixel indices are one-based, whereas SICD is zero-based.
output_meta.ImageData.SCPPixel.Row=uint32(mpdsra_struct.ORP_ROW)-1;
output_meta.ImageData.SCPPixel.Col=uint32(mpdsra_struct.ORP_COLUMN)-1;
% Warning: This assumes PFA.
if isfinite(mpdsra_struct.FOC_X)&&mpdsra_struct.FOC_X~=0&&...
        isfinite(mpdsra_struct.FOC_Y)&&mpdsra_struct.FOC_Y&&...
        isfinite(mpdsra_struct.FOC_Z)&&mpdsra_struct.FOC_Z
    output_meta.PFA.FPN.X=mpdsra_struct.FOC_X;
    output_meta.PFA.FPN.Y=mpdsra_struct.FOC_Y;
    output_meta.PFA.FPN.Z=mpdsra_struct.FOC_Z;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////