function [ output_meta ] = meta2sicd_mensrb( mensrb_struct )
%META2SICD_MENSRB Converts metadata stored in the MENSRB NITF TREs
% into a SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% ARP
lat = latlonnum(mensrb_struct.ACFT_LOC(1:12));
lon = latlonnum(mensrb_struct.ACFT_LOC(13:25));
% Spec says that this value should be MSL, not HAE.  Some sensors populate
% this with HAE anyway.  For these computations (ACFT_ALT and RP_ELV), we
% will attempt to convert from MSL to HAE if possible.  If this value is
% not MSL for a given sensor, then this should be fixed in the
% sensor-specific NITF code.
try % Undulation file may not exist on this system or be on our path
    MSL2HAE=geoid_undulation(lat,lon);
catch
    MSL2HAE=0; % Approximate to ellipsoid
end
ARP_LLA=[lat, lon, (mensrb_struct.ACFT_ALT*FEET_TO_METERS)-MSL2HAE];
ARP=geodetic_to_ecf(ARP_LLA);
output_meta.SCPCOA.ARPPos.X=ARP(1);
output_meta.SCPCOA.ARPPos.Y=ARP(2);
output_meta.SCPCOA.ARPPos.Z=ARP(3);

%% SCP
% MENSRB reference pixel indices are one-based, whereas SICD is zero-based.
output_meta.ImageData.SCPPixel.Row=mensrb_struct.RP_COL-1;
output_meta.ImageData.SCPPixel.Col=mensrb_struct.RP_ROW-1;
lat = latlonnum(mensrb_struct.RP_LOC(1:12));
lon = latlonnum(mensrb_struct.RP_LOC(13:25));
try % Undulation file may not exist on this system or be on our path
    MSL2HAE=geoid_undulation(lat,lon);
catch
    MSL2HAE=0; % Approximate to ellipsoid
end
ORP_LLA=[lat, lon, (mensrb_struct.RP_ELV*FEET_TO_METERS)-MSL2HAE];
output_meta.GeoData.SCP.LLH.Lat=ORP_LLA(1);
output_meta.GeoData.SCP.LLH.Lon=ORP_LLA(2);
output_meta.GeoData.SCP.LLH.HAE=ORP_LLA(3);
ORP=geodetic_to_ecf(ORP_LLA);
output_meta.GeoData.SCP.ECF.X=ORP(1);
output_meta.GeoData.SCP.ECF.Y=ORP(2);
output_meta.GeoData.SCP.ECF.Z=ORP(3);

%% Row/Col unit vectors
row_unit = ned_to_ecf([mensrb_struct.C_R_NC mensrb_struct.C_R_EC mensrb_struct.C_R_DC], ORP, false);
output_meta.Grid.Row.UVectECF=cell2struct(num2cell(row_unit),{'X','Y','Z'});
col_unit = ned_to_ecf([mensrb_struct.C_AZ_NC mensrb_struct.C_AZ_EC mensrb_struct.C_AZ_DC], ORP, false);
output_meta.Grid.Col.UVectECF=cell2struct(num2cell(col_unit),{'X','Y','Z'});

%% SCPCOA
% Not really needed since we can generally derive this from more
% fundamental parameters.
% output_meta.SCPCOA.SlantRange=mensrb_struct.RGCCRP*FEET_TO_METERS;
% output_meta.SCPCOA.GroundRange=mensrb_struct.RGCCRP*mensrb_struct.COSGRZ*FEET_TO_METERS;
% output_meta.SCPCOA.SideOfTrack=mensrb_struct.RLMAP;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////