function add_sar_2kml(k, sicd_meta, varargin)
%ADD_SAR_2KML Adds information about a SAR dataset to an already open KML
%   ADD_SAR_2KML (K, SICD_META, 'PropertyName', PropertyValue, ...)
%   allows property to be set to customize the output.
%
%   SICD_META can be either 1) a SICD structure, 2) a filename of a SAR
%   dataset handled by open_reader() or open_ph_reader, or 3) a reader
%   object as returned by open_reader or open_ph_reader.
%
% TODO: Document property arguments.  For now, see inputParser section of
% code for list and description of possible properties.
%
% Written by: Wade Schwartzkopf, NGA/IDT; Tim Cox, NRL; Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Handle SICD_META argument.  Multiple options here.
if ischar(sicd_meta) && exist(sicd_meta,'file') % Filename was given
    if ~isempty(guess_ph_format(sicd_meta))
        sicd_meta = open_ph_reader(sicd_meta);
    else
        sicd_meta = open_reader(sicd_meta);
    end
    if iscell(sicd_meta)
        for i = 1:numel(sicd_meta)
            add_sar_2kml(k, sicd_meta{i}, varargin{:});
            sicd_meta{i}.close();
        end
        return;
    end
    close_reader_object = true;
else
    close_reader_object = false;
end
if isstruct(sicd_meta) && isfield(sicd_meta,'get_meta') % open_reader object
    reader_obj = sicd_meta;
    sicd_meta = reader_obj.get_meta();
    % If phase history data (metadata is really cphd, rather than sicd)
    if isfield(reader_obj,'read_raw')
        read_fun = reader_obj.read_raw;
    elseif isfield(reader_obj,'read_cphd')
        read_fun = reader_obj.read_cphd;
    end
    if exist('read_fun','var') % Convert CPHD to SICD format
        [ignore, vbmeta] = read_fun('all',[]);
        cphd_meta = sicd_meta;
        sicd_meta = meta2sicd_cphdx(sicd_meta, vbmeta);
    end
end
if ~isfield(sicd_meta,'CollectionInfo')||~isfield(sicd_meta.CollectionInfo,'CoreName')
    sicd_meta.CollectionInfo.CoreName = 'Unknown';
end

% Handle options arguments
p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('srp',true); % Point for spotlight, otherwise line
p.addParamValue('border_color','Red'); % Empty for no border
p.addParamValue('border_thickness',1); % Empty for no border
p.addParamValue('fill_color',[]); % Empty is no fill
p.addParamValue('transparency',1);
p.addParamValue('overlay_max_size',0); % 0 means generate no overlay
p.addParamValue('overlay_decimate','none');
p.addParamValue('overlay_filename',''); % Where to store overlay image.  Path should be relative to kml file location
p.addParamValue('sensor_path',true);
p.addParamValue('sensor_ground_track',false);
p.addParamValue('collection_wedge',true);
p.addParamValue('center_aperture_vector',false);
p.addParamValue('ambiguity_bounds',false); % Currently only in azimuth
p.addParamValue('description_function', @default_SICD_description); % If empty no description is added
p.FunctionName = mfilename;
p.parse(varargin{:});

% Fields that will be used in multiple components
name_str = sicd_meta.CollectionInfo.CoreName;
if isfield(sicd_meta,'Timeline')&&isfield(sicd_meta.Timeline,'CollectStart')
    date_str = datestr(sicd_meta.Timeline.CollectStart,'yyyy-mm-ddTHH:MM:SSZ');
else
    date_str = '';
end
f = k.newFolder(name_str);
% Determine corner lat/longs.  First try sensor model.  Otherwise
% resort to corners in metadata (not guaranteed to be precise).
if isfield(sicd_meta,'ImageData') % True for all complex data, but not phase history
    nx = double(sicd_meta.ImageData.NumCols);
    ny = double(sicd_meta.ImageData.NumRows);
    corner_pix = [1 1 ny ny; 1 nx nx 1];
    corner_latlon = point_slant_to_ground(corner_pix, sicd_meta);
else
    corner_latlon = []; % No corners for phase history
end
if isempty(corner_latlon)
    if isfield(sicd_meta,'GeoData') && ...% Sensor model not available
            isfield(sicd_meta.GeoData,'ImageCorners')&&...
            isfield(sicd_meta.GeoData.ImageCorners,'ICP')
        corner_latlon = zeros(2,4);
        corner_latlon(1,:) = [sicd_meta.GeoData.ImageCorners.ICP.FRFC.Lat, ...
            sicd_meta.GeoData.ImageCorners.ICP.FRLC.Lat, ...
            sicd_meta.GeoData.ImageCorners.ICP.LRLC.Lat, ...
            sicd_meta.GeoData.ImageCorners.ICP.LRFC.Lat];
        corner_latlon(2,:) = [sicd_meta.GeoData.ImageCorners.ICP.FRFC.Lon, ...
            sicd_meta.GeoData.ImageCorners.ICP.FRLC.Lon, ...
            sicd_meta.GeoData.ImageCorners.ICP.LRLC.Lon, ...
            sicd_meta.GeoData.ImageCorners.ICP.LRFC.Lon];
    elseif exist('cphd_meta','var') && isfield(cphd_meta, 'SceneCoordinates') && ...
            isfield(cphd_meta.SceneCoordinates,'ImageAreaCornerPoints') && ...
            isfield(cphd_meta.SceneCoordinates.ImageAreaCornerPoints,'IACP')
        corner_latlon = zeros(2,4);
        corner_latlon(1,:) = [cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(1).Lat, ...
            cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(2).Lat, ...
            cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(3).Lat, ...
            cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(4).Lat];
        corner_latlon(2,:) = [cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(1).Lon, ...
            cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(2).Lon, ...
            cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(3).Lon, ...
            cphd_meta.SceneCoordinates.ImageAreaCornerPoints.IACP(4).Lon];
    end
end
if ~isempty(corner_latlon)
    % Draw image bounding box
    if ~isempty(p.Results.border_color)&&...
        ~isempty(p.Results.border_thickness)&&(p.Results.border_thickness>0)
        kml_args = {'name',name_str,...
                'lineColor',GetKMLColor(p.Results.border_color,1),...
                'lineWidth',p.Results.border_thickness};
        if ~isempty(p.Results.fill_color)
            kml_args = [kml_args {'polyColor',GetKMLColor(p.Results.fill_color, p.Results.transparency)}];
        else
            kml_args = [kml_args {'polyColor','00000000'}];
        end
        if date_str
            kml_args = [kml_args {'timeStamp', date_str}];
        end
        try
            kml_args = [kml_args {'description', p.Results.description_function(sicd_meta)}];
        end
        if size(corner_latlon,1)==3
            f.poly3([corner_latlon(2,:) corner_latlon(2,1)], ...
                [corner_latlon(1,:) corner_latlon(1,1)], ...
                [corner_latlon(3,:) corner_latlon(3,1)], ...
                'altitudeMode','absolute',kml_args{:});
        else
            f.poly([corner_latlon(2,:) corner_latlon(2,1)],...
                [corner_latlon(1,:) corner_latlon(1,1)],...
                'altitudeMode','clampToGround',kml_args{:});
        end
    end
    % Draw overlay image
    if p.Results.overlay_max_size && ~isempty(p.Results.overlay_filename) && ...
            exist('reader_obj','var') && isfield(reader_obj,'read_chip')
        %load chip
        DecSize = p.Results.overlay_max_size;
        DecValue = ceil(max([nx ny])/DecSize);
        chip = double(reader_obj.read_chip([1 nx],[1 ny],[DecValue DecValue],p.Results.overlay_decimate));
        chip_is_real = isreal(chip);
        chip = abs(chip);
        
        % Compute valid data area
        if isfield(sicd_meta.ImageData,'ValidData') && ...
                isfield(sicd_meta.ImageData.ValidData,'Vertex')
            valid_data = poly2mask(...
                round(double([sicd_meta.ImageData.ValidData.Vertex.Row])/DecValue) + 1, ... % ValidData vertices are zero-based
                round(double([sicd_meta.ImageData.ValidData.Vertex.Col])/DecValue) + 1, ... % ValidData vertices are zero-based
                size(chip,1), size(chip,2));
        else
            valid_data = ones(size(chip));
        end

        [cols,rows] = ndgrid(linspace(1,double(sicd_meta.ImageData.NumCols), 5), ...
            linspace(1, double(sicd_meta.ImageData.NumRows), 5)); % Use 5x5 grid
        input_points = [round(rows(:)).'; round(cols(:)).'].';
        latlon = point_slant_to_ground(input_points.', sicd_meta);
        base_points = flipud(latlon(1:2,:)).'; % Want north-south (latitude) in up-down direction
        valid = all(isfinite(base_points),2);
        try
            tform = cp2tform(input_points(valid,:),base_points(valid,:),'polynomial',3);
        catch % Sometimes unable to find the required number of non-collinear points
            try
                tform = cp2tform(input_points(valid,:),base_points(valid,:),'polynomial',2);
            catch
                % Projective ground projection (straight lines remain
                % straight lines).  Not appropriate for most SAR complex
                % imagery, but a last resort.
                % input_points = corner_pix.'; % Really only need 4 points for projective
                % base_points = flipud(corner_latlon(1:2,:)).'; % Want north-south (latitude) in up-down direction
                tform = cp2tform(input_points,base_points,'projective');
            end
        end
        % This will break at the wrapping point for longitude or at the poles.
        north = max(base_points(:,2));
        south = min(base_points(:,2));
        east = max(base_points(:,1));
        west = min(base_points(:,1));
        chip = imtransform(chip, tform, 'size', DecSize.*[1 1], ...
            'FillValues', NaN, 'UData', [1 ny], 'VData', [1 nx], ...
            'XData',[west east],... % Longitude, East right
            'YData',[north south]); % Latitude, North up
        % Probably not necessary to reproject the whole array of the
        % ValidData mask.  Would probably be more efficient to just project
        % ValidData vertices and do poly2mask in the output image space.
        valid_data = imtransform(valid_data, tform, 'size', DecSize.*[1 1], ...
            'FillValues', NaN, 'UData', [1 ny], 'VData', [1 nx], ...
            'XData',[west east],... % Longitude, East right
            'YData',[north south]); % Latitude, North up

        %remap, max decimated typically looks better with "brighter" remap
        if ~chip_is_real
          if strcmpi(p.Results.overlay_decimate,'max')
              chip = brighterremap(chip);
          else
              chip = densityremap(chip);
          end
        else
          chip = uint8(chip);
        end
        
        %write image in path relative to open kml file
        fullfilename = fullfile(fileparts(k.filename),p.Results.overlay_filename);
        imwrite(chip, fullfilename, 'png', 'Alpha', intmax('uint8')*uint8(valid_data));
        
        % Write to KML
        old_includeFiles = k.includeFiles;
        f.overlay(west, east, south, north,...
            'name',name_str,...
            'timeStamp', date_str,...
            'file',p.Results.overlay_filename);
        % Workaround for the fact that the the overlay function in the KML
        % toolbox uses the exact same URL for both the URL reference within
        % the KML and the include files list (which is changed by the
        % overlay function).  We want to use a relative path name within
        % the KML since we might later package all files into a KMZ, but we
        % need the use full URL in the included files list.
        if k.filename
            k.includeFiles = old_includeFiles;
            k.addIncludeFile(fullfilename);
        end
    end
end
if close_reader_object % After overlays we no longer need to read from file
    reader_obj.close();
end

% For now we approximate ARP/GRP paths as straight lines.  Later we can
% add more intermediate steps.
NUM_STEPS = 2;
if isfield(sicd_meta,'ImageFormation') && ...
        all(isfield(sicd_meta.ImageFormation, {'TStartProc','TEndProc'}))
    times = linspace(sicd_meta.ImageFormation.TStartProc,sicd_meta.ImageFormation.TEndProc,2);
elseif exist('vbmeta','var')
    times = linspace(min(vbmeta.TxTime),max(vbmeta.TxTime),NUM_STEPS);
elseif isfield(sicd_meta,'Timeline') && isfield(sicd_meta.Timeline,'CollectDuration')
    times = linspace(0,sicd_meta.Timeline.CollectDuration,2);
else
    times = 0; % Just use center of aperture
end
if isfield(sicd_meta,'Grid') && isfield(sicd_meta.Grid,'TimeCOAPoly')
    coatime = sicd_meta.Grid.TimeCOAPoly(1);
elseif isfield(sicd_meta,'SCPCOA') && isfield(sicd_meta.SCPCOA,'SCPTime')
    coatime = sicd_meta.SCPCOA.SCPTime;
else
    coatime = (times(1) + times(end))/2;
end
    
% Get ARP positions
if isfield(sicd_meta,'Position') && isfield(sicd_meta.Position,'ARPPoly') && ...
        any(times) % Need valid collection times for evaluating polynomial        
    ARPCOAX = polyval(sicd_meta.Position.ARPPoly.X(end:-1:1),coatime);
    ARPCOAY = polyval(sicd_meta.Position.ARPPoly.Y(end:-1:1),coatime);
    ARPCOAZ = polyval(sicd_meta.Position.ARPPoly.Z(end:-1:1),coatime);
    
    ARPPosX = polyval(sicd_meta.Position.ARPPoly.X(end:-1:1),times);
    ARPPosY = polyval(sicd_meta.Position.ARPPoly.Y(end:-1:1),times);
    ARPPosZ = polyval(sicd_meta.Position.ARPPoly.Z(end:-1:1),times);
elseif isfield(sicd_meta,'SCPCOA') && isfield(sicd_meta.SCPCOA,'ARPPos')
    ARPCOAX = sicd_meta.SCPCOA.ARPPos.X;
    ARPCOAY = sicd_meta.SCPCOA.ARPPos.Y;
    ARPCOAZ = sicd_meta.SCPCOA.ARPPos.Z;
    if isfield(sicd_meta.SCPCOA,'ARPVel')
        ARPPosX = sicd_meta.SCPCOA.ARPPos.X + [-1/2 1/2]*(times(end)-times(1))*...
            sicd_meta.SCPCOA.ARPVel.X;
        ARPPosY = sicd_meta.SCPCOA.ARPPos.Y + [-1/2 1/2]*(times(end)-times(1))*...
            sicd_meta.SCPCOA.ARPVel.Y;
        ARPPosZ = sicd_meta.SCPCOA.ARPPos.Z + [-1/2 1/2]*(times(end)-times(1))*...
            sicd_meta.SCPCOA.ARPVel.Z;
    end
end
if exist('ARPPosX','var')
    %convert to LLA
    for i=1:length(ARPPosX)
        pos_lla = ecf_to_geodetic([ARPPosX(i) ARPPosY(i) ARPPosZ(i)]);
        ARPLats(i) = pos_lla(1);
        ARPLons(i) = pos_lla(2);
        ARPAlts(i) = pos_lla(3);
    end
end
if exist('ARPCOAX','var')
    ARPCOA = ecf_to_geodetic([ARPCOAX ARPCOAY ARPCOAZ]);
end

% Get GRP positions
if isfield(sicd_meta,'Position') && isfield(sicd_meta.Position,'GRPPoly') && ...
        any(times) % Need valid collection times for evaluating polynomial
    if isfield(sicd_meta,'CollectionInfo') && ...
            isfield(sicd_meta.CollectionInfo,'RadarMode') && ...
            isfield(sicd_meta.CollectionInfo.RadarMode,'ModeType') && ...
            strcmpi(sicd_meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
        GRPPosX = sicd_meta.Position.GRPPoly.X(1);
        GRPPosY = sicd_meta.Position.GRPPoly.Y(1);
        GRPPosZ = sicd_meta.Position.GRPPoly.Z(1);
        GRPCOAX = GRPPosX;
        GRPCOAY = GRPPosY;
        GRPCOAZ = GRPPosZ;
    else
        GRPPosX = polyval(sicd_meta.Position.GRPPoly.X(end:-1:1),times);
        GRPPosY = polyval(sicd_meta.Position.GRPPoly.Y(end:-1:1),times);
        GRPPosZ = polyval(sicd_meta.Position.GRPPoly.Z(end:-1:1),times);
        GRPCOAX = polyval(sicd_meta.Position.GRPPoly.X(end:-1:1),coatime);
        GRPCOAY = polyval(sicd_meta.Position.GRPPoly.Y(end:-1:1),coatime);
        GRPCOAZ = polyval(sicd_meta.Position.GRPPoly.Z(end:-1:1),coatime);
    end
    %convert to LLA
    for i=1:length(GRPPosX)
        pos_lla = ecf_to_geodetic([GRPPosX(i) GRPPosY(i) GRPPosZ(i)]);
        GRPLats(i) = pos_lla(1);
        GRPLons(i) = pos_lla(2);
        GRPAlts(i) = pos_lla(3);
    end
    GRPCOA = ecf_to_geodetic([GRPCOAX GRPCOAY GRPCOAZ]);
elseif isfield(sicd_meta,'GeoData')&&isfield(sicd_meta.GeoData,'SCP')
    GRPLats = sicd_meta.GeoData.SCP.LLH.Lat;
    GRPLons = sicd_meta.GeoData.SCP.LLH.Lon;
    GRPAlts = sicd_meta.GeoData.SCP.LLH.HAE;
    GRPCOA = [GRPLats GRPLons GRPAlts];
elseif exist('vbmeta','var') && isfield(vbmeta,'SRPPos')
    A = repmat(vbmeta.TxTime, [1 numel(times)]);
    [~, time_ind] = min(abs(A-times));
    [~, coa_time_ind] = min(abs(vbmeta.TxTime-coatime));
    GRPPosX = vbmeta.SRPPos(time_ind,1);
    GRPPosY = vbmeta.SRPPos(time_ind,2);
    GRPPosZ = vbmeta.SRPPos(time_ind,3);
    GRPCOAX = vbmeta.SRPPos(coa_time_ind,1);
    GRPCOAY = vbmeta.SRPPos(coa_time_ind,2);
    GRPCOAZ = vbmeta.SRPPos(coa_time_ind,3);
    %convert to LLA
    for i=1:length(GRPPosX)
        pos_lla = ecf_to_geodetic([GRPPosX(i) GRPPosY(i) GRPPosZ(i)]);
        GRPLats(i) = pos_lla(1);
        GRPLons(i) = pos_lla(2);
        GRPAlts(i) = pos_lla(3);
    end
    GRPCOA = ecf_to_geodetic([GRPCOAX GRPCOAY GRPCOAZ]);
end

if exist('ARPLats','var')
    if p.Results.sensor_path
        a = sprintf('Flight Path for IID: %s',sicd_meta.CollectionInfo.CoreName);
        f.plot3(ARPLons, ARPLats, ARPAlts, 'altitudeMode','absolute',...
            'timeStamp', date_str,...
            'name',a,'description',a,...
            'lineColor',GetKMLColor('Red',1),...
            'lineWidth',p.Results.border_thickness);

        %add start/stop indicators
        f.point(ARPLons(1), ARPLats(1), ARPAlts(1), ...
            'altitudeMode','absolute',...
            'name', 'Start',...
            'description', sprintf('Flight Path Start Point for IID: %s',sicd_meta.CollectionInfo.CoreName),...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Red'));

        f.point(ARPLons(end), ARPLats(end), ARPAlts(end), ...
            'altitudeMode','absolute',...
            'name', 'Stop',...
            'description', sprintf('Flight Path Stop Point for IID: %s',sicd_meta.CollectionInfo.CoreName),...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Red'));
    end

    if p.Results.sensor_ground_track
        a = sprintf('Sensor Ground Track for IID: %s',sicd_meta.CollectionInfo.CoreName);
        f.plot(ARPLons, ARPLats, 'altitudeMode','clampToGround',...
            'timeStamp', date_str,...
            'name',a,'description',a,...
            'lineColor',GetKMLColor('Green',1),...
            'lineWidth',p.Results.border_thickness);

        %add start/stop indicators
        f.point(ARPLons(1), ARPLats(1), 0, ...
            'altitudeMode','clampToGround',...
            'name', 'Start',...
            'description', sprintf('Ground Track Start Point for IID: %s',sicd_meta.CollectionInfo.CoreName),...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Green'));

        f.point(ARPLons(end), ARPLats(end), 0, ...
            'altitudeMode','clampToGround',...
            'name', 'Stop',...
            'description', sprintf('Ground Track Stop Point for IID: %s',sicd_meta.CollectionInfo.CoreName),...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Green'));
    end
end
if exist('GRPLats','var') && p.Results.srp
    if numel(GRPLats)>1
        %add shape
        a = sprintf('Collection GRP Path for IID: %s',sicd_meta.CollectionInfo.CoreName);
        f.plot3(GRPLons, GRPLats, GRPAlts, 'altitudeMode','absolute',...
            'timeStamp', date_str,...
            'name',a,'description',a,...
            'lineColor',GetKMLColor('Cyan',1),...
            'lineWidth',p.Results.border_thickness);

        %add start/stop indicators
        f.point(GRPLons(1), GRPLats(1), GRPAlts(1), ...
            'altitudeMode','absolute',...
            'name', 'Start',...
            'description', sprintf('GRP Path Start Point for IID: %s',sicd_meta.CollectionInfo.CoreName),...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Cyan'));

        f.point(GRPLons(end), GRPLats(end), GRPAlts(end), ...
            'altitudeMode','absolute',...
            'name', 'Stop',...
            'description', sprintf('GRP Path Stop Point for IID: %s',sicd_meta.CollectionInfo.CoreName),...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Cyan'));
    else
        %point mode
        a = sprintf('Collection GRP for IID: %s',sicd_meta.CollectionInfo.CoreName);
        f.point(GRPLons(1), GRPLats(1), GRPAlts(1), ...
            'altitudeMode','absolute',...
            'name', a, 'description', a,...
            'timeStamp', date_str,...
            'iconURL', GetKMLPushpinFile('Cyan'));
    end
end

if exist('ARPLats','var') && exist('GRPLats','var') && p.Results.collection_wedge
    a = sprintf('Collection Wedge for IID: %s',sicd_meta.CollectionInfo.CoreName);
    f.poly3([ARPLons GRPLons(end:-1:1) ARPLons(1)],...
        [ARPLats GRPLats(end:-1:1) ARPLats(1)],...
        [ARPAlts GRPAlts(end:-1:1) ARPAlts(1)],...
        'altitudeMode','absolute',...
        'timeStamp', date_str,...
        'name',a,'description',a,...
        'lineColor',GetKMLColor('White',1),...
        'lineWidth',p.Results.border_thickness,...
        'polyColor',GetKMLColor('White', 0.33));
end
if exist('ARPCOA','var') && exist('GRPCOA','var') && p.Results.center_aperture_vector
    %plot vector from GRP to ARP at center of aperture to
    a = sprintf('Center of Aperture Vector for IID: %s',sicd_meta.CollectionInfo.CoreName);
    f.plot3([GRPCOA(2) ARPCOA(2)], ...
        [GRPCOA(1) ARPCOA(1)], ...
        [GRPCOA(3) ARPCOA(3)], ...
        'altitudeMode','absolute',...
        'timeStamp', date_str,...
        'name',a,'description',a,...
        'lineColor',GetKMLColor('Yellow',1),...
        'lineWidth',p.Results.border_thickness);
end

try % Requires lots of metadata.  Don't want to explicitly check for everything.
    if p.Results.ambiguity_bounds
        %compute where PRF lands for center of scene
        PRFFreq = double(sicd_meta.Timeline.IPP.Set.IPPEnd -  sicd_meta.Timeline.IPP.Set.IPPStart)/...
                      (sicd_meta.Timeline.IPP.Set.TEnd - sicd_meta.Timeline.IPP.Set.TStart);
        dca = sicd_meta.SCPCOA.DopplerConeAng;
        GRP = [sicd_meta.GeoData.SCP.ECF.X sicd_meta.GeoData.SCP.ECF.Y sicd_meta.GeoData.SCP.ECF.Z];
        ARV = [sicd_meta.SCPCOA.ARPVel.X; sicd_meta.SCPCOA.ARPVel.Y; sicd_meta.SCPCOA.ARPVel.Z];
        ARP = [sicd_meta.SCPCOA.ARPPos.X sicd_meta.SCPCOA.ARPPos.Y sicd_meta.SCPCOA.ARPPos.Z];
        R = norm(abs(GRP-ARP));    
        V = norm(ARV);
        rov = (R/V);
        thetaDot = sind(dca)/rov;
        vfrqC = sicd_meta.Grid.Row.KCtr*SPEED_OF_LIGHT/2;

        vd = (SPEED_OF_LIGHT*PRFFreq)/(2*thetaDot*vfrqC);
        PixelDist = round(vd/(2*sicd_meta.Grid.Col.SS));

        %now form left line in image space and convert to LLA
        CenterX = round(sicd_meta.ImageData.NumCols/2);
        points = [1 CenterX-PixelDist; sicd_meta.ImageData.NumRows CenterX-PixelDist];
        %convert to lla
        pos = point_slant_to_ground(points', sicd_meta);

        %write shapefile
        a = sprintf('Left PRF Line for IID: %s',sicd_meta.CollectionInfo.CoreName);
        f.plot3([pos(2,1) pos(2,2)], ...
            [pos(1,1) pos(1,2)], ...
            [pos(3,1) pos(3,2)], ...
            'altitudeMode','absolute',...
            'timeStamp', date_str,...
            'name',a,'description',a,...
            'lineColor',GetKMLColor('Green',1),...
            'lineWidth',p.Results.border_thickness);

        %now form right line in image space and convert to LLA
        CenterX = round(sicd_meta.ImageData.NumCols/2);
        points = [1 CenterX+PixelDist; sicd_meta.ImageData.NumRows CenterX+PixelDist];
        %convert to lla
        pos = point_slant_to_ground(points', sicd_meta);

        %write shapefile
        a = sprintf('Right PRF Line for IID: %s',sicd_meta.CollectionInfo.CoreName);
        f.plot3([pos(2,1) pos(2,2)], ...
            [pos(1,1) pos(1,2)], ...
            [pos(3,1) pos(3,2)], ...
            'altitudeMode','absolute',...
            'timeStamp', date_str,...
            'name',a,'description',a,...
            'lineColor',GetKMLColor('Green',1),...
            'lineWidth',p.Results.border_thickness);
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////