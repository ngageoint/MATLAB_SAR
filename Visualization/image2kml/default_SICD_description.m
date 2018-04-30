function str = default_SICD_description(sicdmeta, Units)
%DEFAULT_SICD_DESCRIPTION Converts a SICD structure to a string
%  DEFAULT_SICD_DESCRIPTION(sicdmeta) takes
%  a SICD-compatible structure and generates a Google Earth compatible KML
%  Placemark entity description (HTML).
%
% Note: No this isn't the most beautiful of HTML...  We (i.e. you) should
% probably provide a description table style and link to it and include the
% specific information in a more usable format for your particular use.
% This is meant to be a stub...
%
% Author: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('Units','var')
    Units = 'Metric';
end

str = sprintf('<table width="90%%" border="1">\n');
if isfield(sicdmeta,'CollectionInfo')
    if isfield(sicdmeta.CollectionInfo,'Classification')
        str = [str sprintf('<tr><th>Classification:</th>              <td>%s</td></tr>\n',sicdmeta.CollectionInfo.Classification)];
    end
    if isfield(sicdmeta.CollectionInfo,'CoreName')
        str = [str sprintf('<tr><th>Core Name:</th>                   <td>%s</td></tr>\n',sicdmeta.CollectionInfo.CoreName)];
    end
    if isfield(sicdmeta.CollectionInfo,'CollectorName')
        str = [str sprintf('<tr><th>Collector Name:</th>              <td>%s</td></tr>\n',num2str(sicdmeta.CollectionInfo.CollectorName))];
    end
    if isfield(sicdmeta.CollectionInfo,'RadarMode')
        str = [str '<tr><th>Mode:</th>                        <td>'];
        if isfield(sicdmeta.CollectionInfo.RadarMode,'ModeID')
            str = [str sicdmeta.CollectionInfo.RadarMode.ModeID];
        end
        if all(isfield(sicdmeta.CollectionInfo.RadarMode,{'ModeID','ModeType'}))
            str = [str ' / '];
        end
        if isfield(sicdmeta.CollectionInfo.RadarMode,'ModeType')
            str = [str sicdmeta.CollectionInfo.RadarMode.ModeType];
        end
        str = [str sprintf('</td></tr>\n')];
    end
end
if isfield(sicdmeta,'ImageFormation')&&isfield(sicdmeta.ImageFormation,'TxRcvPolarizationProc')&&...
        ischar(sicdmeta.ImageFormation.TxRcvPolarizationProc)
    PolString = sicdmeta.ImageFormation.TxRcvPolarizationProc;
    str = [str sprintf('<tr><th>Polarization:</th>                <td>%c%c</td></tr>\n',PolString(1),PolString(3))];
end
if isfield(sicdmeta,'Timeline')&&isfield(sicdmeta.Timeline,'CollectStart')
    ZuluTime = [datestr(sicdmeta.Timeline.CollectStart,'dd-mmm-yyyy HH:MM') ' Z'];
    if isfield(sicdmeta.GeoData,'SCP')
        days_past_zulu=sicdmeta.GeoData.SCP.LLH.Lon/360;
        LSTTime = [datestr(sicdmeta.Timeline.CollectStart+days_past_zulu,'HH:MM') ' LS']; % Local solar time
    else
        LSTTime = '';
    end
    str = [str sprintf('<tr><th>Collection Date:</th>             <td>%s / % s</td></tr>\n', ZuluTime, LSTTime)];
end
if isfield(sicdmeta,'ImageData')&&all(isfield(sicdmeta.ImageData,{'NumRows','NumCols'}))
    str = [str sprintf('<tr><th>Image Size (rows x cols):</th>    <td>%d x %d (px)</td></tr>\n',sicdmeta.ImageData.NumRows,sicdmeta.ImageData.NumCols)];
end
if isfield(sicdmeta,'Grid')&&all(isfield(sicdmeta.Grid,{'Row','Col'}))
    if isfield(sicdmeta.Grid.Row,'ImpRespWid')&&isfield(sicdmeta.Grid.Col,'ImpRespWid')        
        if strcmp(Units,'Metric')
            str = [str sprintf('<tr><th>Resolution (rng x az):</th>    <td>%.2f x %.2f (m)</td></tr>\n',sicdmeta.Grid.Row.ImpRespWid,sicdmeta.Grid.Col.ImpRespWid)];
        else
            if (sicdmeta.Grid.Row.ImpRespWid/FEET_TO_METERS < 1 && sicdmeta.Grid.Col.ImpRespWid/FEET_TO_METERS < 1)
                str = [str sprintf('<tr><th>Resolution (rng x az):</th>    <td>%.2f x %.2f (inches)</td></tr>\n',12*sicdmeta.Grid.Row.ImpRespWid/FEET_TO_METERS,12*sicdmeta.Grid.Col.ImpRespWid/FEET_TO_METERS)];
            else
                str = [str sprintf('<tr><th>Resolution (rng x az):</th>    <td>%.2f x %.2f (ft)</td></tr>\n',sicdmeta.Grid.Row.ImpRespWid/FEET_TO_METERS,sicdmeta.Grid.Col.ImpRespWid/FEET_TO_METERS)];
            end
        end
    end
    if isfield(sicdmeta.Grid.Row,'SS')&&isfield(sicdmeta.Grid.Col,'SS')
        if strcmp(Units,'Metric')
            str = [str sprintf('<tr><th>Sample spacing (rng x az):</th>    <td>%.2f x %.2f (m)</td></tr>\n',sicdmeta.Grid.Row.SS,sicdmeta.Grid.Col.SS)];
        else
            if (sicdmeta.Grid.Row.ImpRespWid/FEET_TO_METERS < 1 && sicdmeta.Grid.Col.ImpRespWid/FEET_TO_METERS < 1)
                str = [str sprintf('<tr><th>Sample spacing (rng x az):</th>    <td>%.2f x %.2f (inches)</td></tr>\n',12*sicdmeta.Grid.Row.SS/FEET_TO_METERS,12*sicdmeta.Grid.Col.SS/FEET_TO_METERS)];
            else
                str = [str sprintf('<tr><th>Sample spacing (rng x az):</th>    <td>%.2f x %.2f (ft)</td></tr>\n',sicdmeta.Grid.Row.SS/FEET_TO_METERS,sicdmeta.Grid.Col.SS/FEET_TO_METERS)];
            end
        end
    end
end
if isfield(sicdmeta,'SCPCOA')
    if isfield(sicdmeta.SCPCOA,'SideOfTrack')
        str = [str sprintf('<tr><th>Look Side:</th>                   <td>%s</td></tr>\n',sicdmeta.SCPCOA.SideOfTrack)];
    end
    if isfield(sicdmeta.SCPCOA,'GrazeAng')
        str = [str sprintf('<tr><th>Grazing Angle:</th>               <td>%.1f (deg)</td></tr>\n',sicdmeta.SCPCOA.GrazeAng)];
    end
    if isfield(sicdmeta.GeoData,'SCP')&&isfield(sicdmeta.SCPCOA,'ARPPos')
        GRP = [sicdmeta.GeoData.SCP.ECF.X sicdmeta.GeoData.SCP.ECF.Y sicdmeta.GeoData.SCP.ECF.Z];
        ARP = [sicdmeta.SCPCOA.ARPPos.X sicdmeta.SCPCOA.ARPPos.Y sicdmeta.SCPCOA.ARPPos.Z];
        gpn = wgs_84_norm(GRP); % Ground plane normal, based on wgs_84 ellipsoid
        range = ARP-GRP; % Range vector
        north_ground = [ 0 0 1 ]-([ 0 0 1 ]*gpn)*gpn.'; % Project north onto ground plane
        range_ground=range-(range*gpn)*gpn.'; % Project range vector onto ground plane
        azimuth = atan2( dot( cross( range_ground, north_ground ), gpn ), dot( north_ground, range_ground ) );
        azimuth = mod(azimuth*180/pi,360); % Degrees in [0,360] range
        str = [str sprintf('<tr><th>Azimuth Angle:</th>               <td>%.1f (deg)</td></tr>\n',azimuth)];
    end
    if isfield(sicdmeta.SCPCOA,'DopplerConeAng')
        str = [str sprintf('<tr><th>Doppler Cone Angle:</th>          <td>%.1f (deg)</td></tr>\n',sicdmeta.SCPCOA.DopplerConeAng)];
    end
end
str = [str '</table>'];

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////