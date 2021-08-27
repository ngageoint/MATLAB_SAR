function grid_params = bp_parse_grid_params( per_pulse_metadata, num_samples_used, varargin )
%BP_PARSE_GRID_PARAMS Compute grid for backprojection
% BP_PARSE_GRID_PARAMS(PER_PULSE_METADATA, NUM_SAMPLES_USED, ...
%    'PropertyName', PropertyValue, ...)
%
%       Property name     Description
%       center            Center of image to form in ECF coordinates.
%                         Default is scene reference point.
%       grid_type         'slant', 'ground', or 'DEM'.  Default is ground.
%       image_size_meters Scene size for image formation ([range_size
%                         azimuth_size]) in meters.  Default is total size
%                         supported by collect.
%       image_size_pixels Size of the image in pixels (only used for images
%                         formed to a DEM)
%       sample_rate       Samples per IPR.  Default is 1.5.
%       grid              The grid of image points (ECEF) onto which the
%                         image should be formed.  The format of the grid
%                         parameters is identical to the (MxNx3) grid
%                         returned from this routine if no grid is
%                         supplied.
%
% Output is a structure with the following fields:
%       center            Same as input property
%       grid_type         Same as input property
%       grid              Same as input property
%       sample_rate       Same as input property
%       range_vect_hat    Unit vector point from center of aperture to
%                         image center
%       image_size_meters Extent of the image in meters
%       image_size_pixels Size of the image in pixels
%       row_coords        Vector describing position of each row.  For
%                         planar grids, this is the distance from the
%                         center of the image in meters.  For DEMs, this is
%                         the latitude position in degrees of each row.
%       col_coords        Vector describing position of each column.
%       row_unit_vector   For planar grids, the unit vector in direction of
%                         increasing row
%       col_unit_vector   For planar grids, the unit vector in direction of
%                         increasing column
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Calculate default image extent and supported resolution (in slant plane)
pulse_bandwidth = per_pulse_metadata.SCSS * (num_samples_used-1);
[max_resolution, max_extent] = pulse_info_to_resolution_extent(...
    per_pulse_metadata.TxPos - per_pulse_metadata.SRPPos, ... % Line-of-sight vector between ARP and ORP
    per_pulse_metadata.SC0 + (pulse_bandwidth/2), per_pulse_metadata.SCSS, ...
    pulse_bandwidth);

% Parse input parameters
p1 = inputParser;
p1.KeepUnmatched=true;
p1.addParamValue('grid_type',  'ground',                     @(x) any(strcmp(x,{'slant','ground','DEM'})));
p1.addParamValue('grid',       []); % Single grid that fits in memory.  Perhaps in the future, this could also allow for a file input.
p1.addParamValue('center',     mean(per_pulse_metadata.SRPPos), @(x) isempty(x) | (isvector(x)&&(length(x)==3))); % Centerpoint of formed image (ECEF, meters)
p1.addParamValue('image_size_meters', max_extent,            @(x) (isvector(x)&&(length(x)==2))); % Scene extent [x y] (m)
p1.addParamValue('image_size_pixels', [],                    @(x) (isvector(x)&&(length(x)==2))); % Scene extent [x y] (pixels)
p1.addParamValue('sample_rate',1.5,                          @(x) (isvector(x)&&(length(x)<=2)));
p1.FunctionName = mfilename;
p1.parse(varargin{:});
grid_params = p1.Results;

% Some geometry info
range_vect = grid_params.center - ... % Range vector points from center of aperture to image center
    ((per_pulse_metadata.TxPos(ceil(length(per_pulse_metadata.TxPos)/2),:) + ... % Use middle pulse
    per_pulse_metadata.RcvPos(ceil(length(per_pulse_metadata.RcvPos)/2),:))/2); % Bistatic definition, average of tx and rcv
% range_vect = grid_params.center - ... % Range vector points from center of aperture to image center
%     per_pulse_metadata.TxPos(ceil(length(per_pulse_metadata.TxPos)/2),:); % Use middle pulse
grid_params.range_vect_hat = range_vect/norm(range_vect); % Unit vector
vel_vector = per_pulse_metadata.TxPos(end,:)-per_pulse_metadata.TxPos(1,:); % Velocity vector
spn=cross(grid_params.range_vect_hat,vel_vector); % Slant plane normal
spn=spn/norm(spn); % Unit vector
% Slant plane normal to point "up" (away from origin).  This is affected by left/right looking.
spn = spn * sign(grid_params.center*cross(grid_params.range_vect_hat,vel_vector).'); % Left/right fix

% Determine final image size.  IMAGE_SIZE_PIXELS input parameter is ignored
% in all cases except for DEM.  We only allow it for the DEM case, since we
% don't currently have a robust automated way to determine it there.
if ~isempty(grid_params.grid) % Arbitrary grid was passed in.  Nothing to compute.
    grid_params.grid_type = 'ignore';
    grid_params.center = squeeze(mean(mean(grid_params.grid,1),2)); % Centroid of points in grid
    % grid_params.image_size_meters = ???
    grid_params.image_size_pixels = [size(grid_params.grid,1), size(grid_params.grid,2)];
    grid_params = rmfield(grid_params,'sample_rate'); % Not valid
    return; % Nothing left to do for explicitly passed grid
elseif strcmpi(grid_params.grid_type,'DEM')
    if isempty(grid_params.image_size_pixels) % Default
        % Assures that at least sample rate is achieved in all directions.
        grid_params.image_size_pixels = ceil(grid_params.image_size_meters.*...
            max(grid_params.sample_rate)./min(max_resolution));
    end
else % Planar image (slant or ground)
    % Adjust distances from slant to ground if necessary
    if strcmpi(grid_params.grid_type,'ground')
        gpn = wgs_84_norm( grid_params.center ); % Ground plane normal, based on (inflated) wgs_84 ellipsoid
        graze = asin(gpn.'*(-grid_params.range_vect_hat.')); % radians
        twist = asin(cross(gpn.',spn)*(-grid_params.range_vect_hat.')/...
            norm(cross(-grid_params.range_vect_hat.',gpn.'))); % radians
        if any(strcmp(p1.UsingDefaults,'image_size_meters'))
            grid_params.image_size_meters(1) = grid_params.image_size_meters(1)/cos(graze);
            grid_params.image_size_meters(2) = grid_params.image_size_meters(2)/cos(twist);
        end        
        max_resolution(1) = max_resolution(1)/cos(graze);
        max_resolution(2) = max_resolution(2)/cos(twist);
    end
    grid_params.image_size_pixels = ceil(grid_params.image_size_meters.*...
        grid_params.sample_rate./max_resolution);
end

% Compute grid onto which the image will be formed
switch grid_params.grid_type
    % For planar images (ground and slant) compute the unit vectors in the
    % direction of increasing row/column direction
    case 'ground'
        row_vector = grid_params.range_vect_hat-(grid_params.range_vect_hat*gpn)*gpn.'; % Project range vector to ground
        grid_params.row_unit_vector = row_vector/norm(row_vector); % Unit vector
        col_vector = -cross(grid_params.row_unit_vector,gpn); % Vector in direction of increasing column (ECEF)
        grid_params.col_unit_vector = col_vector/norm(col_vector); % Unit vector
    case 'slant'
        % Slant plane really has no meaning in backprojection.  Any
        % arbitrary grid can be defined.  We define it here merely for
        % comparison against PFA.
        grid_params.row_unit_vector = grid_params.range_vect_hat; % Unit vector in direction of increasing row (ECEF)
        grid_params.col_unit_vector = -cross(grid_params.row_unit_vector,spn); % Unit vector in direction of increasing column (ECEF)
    case 'DEM'
        % We'll also get the DEM data for the image.  We assume that the
        % entire image fits within .2 degrees of lat and long and that the
        % entire DEM can fit into memory (reasonable for SRTM, but probably
        % not LIDAR).
        center_lla  = ecf_to_geodetic(grid_params.center);
        [lats,lons,elevations] = get_DEM_heights_region([center_lla(1)-0.1,center_lla(2)-0.1], ...
            [center_lla(1)+0.1,center_lla(2)+0.1]);
        [lats, lons] = meshgrid(lats, lons);
        % Return DEM interpolation function, rather than the grid itself,
        % since for image formation we will need values not lying exactly
        % on the DEM points.
        grid_params.F_1 = TriScatteredInterp(lats(:), lons(:), elevations(:));
        
        % Convert scene extent from meters to arc-degrees
        e2 = 6.6943799901377997e-3;  % eccentricity squared of Earth (WGS 84 value)
        a = 6378137.0;               % semimajor radius of the Earth (WGS 84 value)
        b = a .* sqrt( 1.0 - e2);    % semiminor radius of the Earth
        lat = pi/180*center_lla(1);
        M = (a*b)^2/((a*cos(lat))^2+(b*sin(lat))^2)^(3/2);
        N = a^2/sqrt((a*cos(lat))^2+(b*sin(lat))^2);
        arcradius_lat = pi/180*M;          % length of arc-degree in latitude direction (m)
        arcradius_lon = pi/180*cos(lat)*N; % length of arc-degree in longitude direction (m)
        scene_extent_lat = grid_params.image_size_meters(1)/arcradius_lat;
        scene_extent_lon = grid_params.image_size_meters(2)/arcradius_lon;
end
% Compute row/column coordinates
switch grid_params.grid_type
    % For a plane, coordinates are distances (in meters) from the center
    % of image for each row/column
    case {'ground','slant'}
        grid_params.row_coords = grid_params.image_size_meters(1)*...
            linspace(-1,1,grid_params.image_size_pixels(1))/2;
        grid_params.col_coords = grid_params.image_size_meters(2)*...
            linspace(-1,1,grid_params.image_size_pixels(2))/2;
    % For a DEM, coordinates of each row/column are in lat/lon degrees;
    % that is, angularly, rather than along a plane.
    case 'DEM'
        grid_params.row_coords = center_lla(1) + (scene_extent_lat*...
            linspace(-1,1,grid_params.image_size_pixels(1))/2);
        grid_params.col_coords = center_lla(2) + (scene_extent_lon*...
            linspace(-1,1,grid_params.image_size_pixels(2))/2);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////