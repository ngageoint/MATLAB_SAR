function [ sicd_meta ] = meta2sicd_palsar2img( native_meta )
%META2SICD_PALSAR2IMG Converts ALOS PALSAR 2 image file metadata into a SICD-style metadata structure
%
% Takes as input a metadata structure from read_ceos_img_meta.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% ImageData
sicd_meta.ImageData.NumRows=uint32(native_meta.num_pixels);
sicd_meta.ImageData.NumCols=uint32(native_meta.num_lines);
sicd_meta.ImageData.FullImage=sicd_meta.ImageData; % Assume full image is given
if strncmp(native_meta.sar_datatype_code, 'C*8', 3)  % Should be true for all level 1.1
    sicd_meta.ImageData.PixelType='RE32F_IM32F';
end
sicd_meta.ImageData.FirstRow=uint32(0);
sicd_meta.ImageData.FirstCol=uint32(0);

%% GeoData
% Just seeding with a rough value.  GeodData.SCP will be computed precisely later.
sicd_meta.GeoData.SCP.LLH.Lat = mean([native_meta.signal.lat_center])/1e6;
sicd_meta.GeoData.SCP.LLH.Lon = mean([native_meta.signal.lon_center])/1e6;
% Average terrain height above elliposid at scene center is not provided in level 1.1 data
sicd_meta.GeoData.SCP.LLH.HAE = 0;  % Perhaps we should pull a real value from an external DEM source?
ecf=geodetic_to_ecf([sicd_meta.GeoData.SCP.LLH.Lat sicd_meta.GeoData.SCP.LLH.Lon sicd_meta.GeoData.SCP.LLH.HAE]);
sicd_meta.GeoData.SCP.ECF.X=ecf(1);
sicd_meta.GeoData.SCP.ECF.Y=ecf(2);
sicd_meta.GeoData.SCP.ECF.Z=ecf(3);

%% Timeline
% Time of first line, not necessarilly first transmitted pulse.
sicd_meta.Timeline.CollectStart = datenum(native_meta.signal(1).year, ...
    0, native_meta.signal(1).day, 0, 0, native_meta.signal(1).usec/1e6);
% Not really collect duration, but rather range of line times
% CollectDuration is required in SICD so we must make up something.  We
% could add the dwell time to this, and this would probably a bit more
% accurate.
sicd_meta.Timeline.CollectDuration = (native_meta.signal(end).usec - ...
    native_meta.signal(1).usec) * 1e-6;
sicd_meta.Timeline.IPP.Set.TStart = 0;
sicd_meta.Timeline.IPP.Set.TEnd = sicd_meta.Timeline.CollectDuration;
sicd_meta.Timeline.IPP.Set.IPPStart=uint32(0);
prf = native_meta.signal(1).prf/1000;
sicd_meta.Timeline.IPP.Set.IPPPoly=[0; prf];
sicd_meta.Timeline.IPP.Set.IPPEnd = ...
    uint32(polyval(sicd_meta.Timeline.IPP.Set.IPPPoly(end:-1:1),...
    double(sicd_meta.Timeline.IPP.Set.TEnd)));

%% RadarCollection
if native_meta.signal(1).tx_pol(1)==0
    tx_pol = 'H';
else
    tx_pol = 'V';
end
if native_meta.signal(1).rcv_pol(1)==0
    rcv_pol = 'H';
else
    rcv_pol = 'V';
end
% Might be others as well, but this is the TX polarization for this channel
sicd_meta.RadarCollection.TxPolarization = tx_pol;

%% ImageFormation
sicd_meta.ImageFormation.TxRcvPolarizationProc = [tx_pol ':' rcv_pol];
% We use these times only because we have nothing better.
sicd_meta.ImageFormation.TStartProc = 0;
sicd_meta.ImageFormation.TEndProc = sicd_meta.Timeline.CollectDuration;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////