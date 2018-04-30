function [ output_meta ] = meta2sicd_gff( input_meta )
%META2SICD_GFF Converts metadata from read_gff_meta into SICD-style structure
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Even GFFs with datatype of INT are returned as 32-bit IQ, since INT types
% are store phase-magnitude in GFF.
output_meta.ImageData.PixelType='RE32F_IM32F';
output_meta.ImageData.NumCols=uint32(input_meta.AzCnt);
output_meta.ImageData.NumRows=uint32(input_meta.RgCnt);
output_meta.ImageData.FullImage=output_meta.ImageData; % Assume full image
output_meta.ImageData.FirstRow=uint32(0); output_meta.ImageData.FirstCol=uint32(0);
output_meta.ImageData.SCPPixel.Row=output_meta.ImageData.NumRows/2;
output_meta.ImageData.SCPPixel.Col=output_meta.ImageData.NumCols/2;

output_meta.GeoData.EarthModel='WGS_84'; % Constant for all SICD
output_meta.GeoData.SCP.LLH.Lat=input_meta.SRPLat;
output_meta.GeoData.SCP.LLH.Lon=input_meta.SRPLong;
output_meta.GeoData.SCP.LLH.HAE=input_meta.SRPAlt;
srp_ecf = geodetic_to_ecf([input_meta.SRPLat input_meta.SRPLong input_meta.SRPAlt]);
output_meta.GeoData.SCP.ECF.X=srp_ecf(1);
output_meta.GeoData.SCP.ECF.Y=srp_ecf(2);
output_meta.GeoData.SCP.ECF.Z=srp_ecf(3);

output_meta.Grid.Row.SS=input_meta.RgPixelSz;
output_meta.Grid.Row.ImpRespWid=input_meta.RgResolution;
% No info in metadata on zeropad, so we assume Nyquist sampling.  Not really true though.
% If data was uniformly weighted, we could compute ImpRespBW from
% ImpRespWid, but no information is given on weighting, so we either have
% to guess or determine from data.
output_meta.Grid.Row.ImpRespBW=1/output_meta.Grid.Row.SS;
output_meta.Grid.Row.DeltaK1=output_meta.Grid.Row.ImpRespBW/2;
output_meta.Grid.Row.DeltaK2=-output_meta.Grid.Row.DeltaK1;
output_meta.Grid.Row.DeltaKCOAPoly = 0.5/output_meta.Grid.Row.SS; % Data essentially has FFTSHIFT applied
output_meta.Grid.Row.Sgn=-1;
output_meta.Grid.Col.SS=input_meta.AzPixelSz;
output_meta.Grid.Col.ImpRespWid=input_meta.AzResolution;
% No info in metadata on zeropad, so we assume Nyquist sampling.  Not really true though.
output_meta.Grid.Col.ImpRespBW=1/output_meta.Grid.Col.SS;
output_meta.Grid.Col.DeltaK1=output_meta.Grid.Col.ImpRespBW/2;
output_meta.Grid.Col.DeltaK2=-output_meta.Grid.Col.DeltaK1;
output_meta.Grid.Col.DeltaKCOAPoly = 0.5/output_meta.Grid.Col.SS; % Data essentially has FFTSHIFT applied
output_meta.Grid.Col.Sgn=-1;

output_meta.SCPCOA.GrazeAng=input_meta.GrzAngle;
output_meta.SCPCOA.IncidenceAng=90-output_meta.SCPCOA.GrazeAng;
if input_meta.Squint<0
    output_meta.SCPCOA.SideOfTrack='L';
else
    output_meta.SCPCOA.SideOfTrack='R';
end

% Could add center frequency input_meta.NominalCenterFreq
%    Would need to get bandwidth from resolution
% Could get SCPCOA.ARPPos from converting APCLat, APCLong, and APCAlt
% APC, velocity (V), and SRP in ECEF are fields in GFF, but are left blank
%    in the examples known about.

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////