function undulation = geoid_undulation(lat, lon, undulationFilename)
% GEOID_UNDULATION Computes the difference between the WGS84
% ellipsoid and the EGM2008 Global Gravitational geoid.  This
% "undulation" is computed by interpolating between point values
% known on a grid.  The point (undulation) values are read
% from a file of format described in "README_WGS84_2.pdf".
%
% This function opens and closes the specified file to extract the required
% undulation data surrounding the specified point.  It is probably not
% efficient enough to correct large amounts of data.
%
% Reference:
%   http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
%
% Inputs:
%       lat:            A floating point decimal degree value for the
%                       latitude of the point whos undulation is to be
%                       determined.
%       lon:            A floating point decimal degree value for the
%                       longitude of the point whos undulation is to be
%                       determined.
%
% Allowed properties:
%       Property name         Description
%       -----------------------------------------------------------------
%       undulationFilename The name of the geoid undulation file that
%                          contains elevation offsets.  If a file name is
%                          provided it will be used; otherwise a default
%                          file name will be searched for on MATLAB path.
%
% Written by: Tom Krauss, NGA/IDT
%             Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

defaultFilenames = {...
    'Und_min1x1_egm2008_isw=82_WGS84_TideFree_SE', ...
    'Und_min1x1_egm2008_isw=82_WGS84_TideFree', ...
    'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE', ...
    'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree'};

if ~exist('undulationFilename', 'var') || isempty(undulationFilename)
    % Look to see if one of the default filenames is on the MATLAB path
    for i = 1:numel(defaultFilenames)
        % Note: we need to append the '.' to the file name because "which" uses
        % that to identify the string as a file (rather than a function).  If
        % the file name already has a '.' it still works, but if it's missing
        % one "which" will fail.
        undulationFilename = which([defaultFilenames{i} '.']);
        if ~isempty(undulationFilename), break, end
    end
    if isempty(undulationFilename)
        error('GEOID_UNDULATION:UNDULATION_FILE_NOT_FOUND',...
            'Undulation file not found.  Pass filename explicitly or place a file with the default filename on the MATLAB path.');
    end
end

fileStats = dir(undulationFilename);
b = fileStats.bytes;
% Please excuse the hard-coded numbers here.  This is solving the
% inverse function mapping grid size to number of bytes.  That is,
% solve b = (180/(x/60) + 1)(360/(x/60)) for x (here x=grid).  The
% hard-coded numbers are a result of the 180, 360, 60, etc.
grid = (-172800 - sqrt(172800^2-4*(8-b)*933120000))/(2*(8-b));

% The file is written using Fortran records so the record length
% includes a 4-byte "leader" and a 4-byte "trailer" delimiting the
% records.  Each undulation value is a 4-byte float.
delta_lat          = grid/60;
delta_lon          = grid/60;
num_records        = 180/delta_lat;          % -90 to +90 in delta_lat divisions
num_lon_values     = 360/delta_lon;          % 0 to 360 in delta_lon divisions
num_record_entries = num_lon_values+2;       % num. divisions plus record delimiters
record_length      = 4*(num_record_entries); % 4-byte float values

% The grid requires 0 <= longitude < 360 so we'll "unwrap" negative
% values.
lon = mod(lon,360);

% Read range of values that span requested points
% This functions grabs undulations values as a single read from a
% rectangular subarray. If requested points are widely dispersed, it may be
% more efficient to do multiple smaller reads from across the file.
% All the min/max's in the index computations below come from the fact that
% you must grab at least 2 lat/lons for interp2 to work, even if the point
% for evaluation falls exactly on a sample.
lat_min_ind = min(floor(1+(90-max(lat(:)))/delta_lat),num_records-1);
lat_max_ind = max(ceil(1+(90-min(lat(:)))/delta_lat),lat_min_ind+1);
lon_min_ind = floor(1+min(lon(:))/delta_lon);
lon_max_ind = min(max(ceil(1+max(lon(:))/delta_lon),lon_min_ind+1),num_lon_values);
latsize=lat_max_ind-lat_min_ind+1;
lonsize=lon_max_ind-lon_min_ind+1;

% Read relevant values from file
fid = fopen(undulationFilename, 'r', 'l');
% Detect endianness
if fread(fid,1,'*uint32') ~= num_lon_values*4
    fclose(fid);
    fid = fopen(undulationFilename, 'r', 'b');
    if fread(fid,1,'*uint32') ~= num_lon_values*4
        error('GEOID_UNDULATION:INVALID_FILE',...
            'Undulation file not in expected format.');
    end
end
fseek(fid,4+... % Beginning of data
    ((lat_min_ind-1)*record_length)+... % Skip to first row of interest
    ((lon_min_ind-1)*4),'bof'); % Skip to first column of interest
vals=fread(fid,[lonsize,latsize],[num2str(lonsize) '*float32'],...
    8+(num_lon_values-lonsize)*4)';
if max(lon(:)) > (360-delta_lon) % Replicate 0 for 360
    fseek(fid,4+((lat_min_ind-1)*record_length),'bof'); % First longitude
    vals=[vals, fread(fid,[1 latsize],'1*float32',...
        8+(num_lon_values-lonsize)*4)'];
    lon_max_ind = lon_max_ind + 1;  % For later interpolation
end
fclose(fid);

% We'll use Matlab's built-in interpolator so we'll need to build arrays to
% represent the latitude and longitude "grid" of the points.  The "Xs" and
% "Ys" arrays give the latitudes and longitudes of the requested points
% while "vals" gives the undulation at those points.
[Xs,Ys] = meshgrid(delta_lon*((lon_min_ind:lon_max_ind)-1),...
    90 - delta_lat*((lat_min_ind:lat_max_ind)-1));
undulation = interp2(Xs, Ys, vals, lon, lat, 'linear');
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////