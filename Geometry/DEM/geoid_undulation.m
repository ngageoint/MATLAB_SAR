function undulation = geoid_undulation(lat, lon, varargin)
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
%       PropertyName:   The PropertyName and PropertyValue are key-value
%       PropertyValue:  pairs used throughout the tool kit...
%
% Allowed properties:
%       Property name         Description
%       -----------------------------------------------------------------
%       undulationFilename The name of the geoid undulation file that
%                          contains elevation offsets.  If a file name is
%                          provided it will be used otherwise a default
%                          file will be used.
%       useHiRes           If a file name is not provided via the
%                          undulationFilename parameter a default will be
%                          used.  The default is the "low resolution"
%                          2.5x2.5 minute grid.  If the useHiRes parameter
%                          is set to true, however, the higher-resolution
%                          1x1 minute grid will be used.
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

defaultLowResUundulationFilename = 'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE';
defaultHiResUundulationFilename  = 'Und_min1x1_egm2008_isw=82_WGS84_TideFree_SE';


p = inputParser;
p.addParamValue('undulationFilename', defaultLowResUundulationFilename);
p.addParamValue('useHiRes', false);
p.parse(varargin{:});

filename = p.Results.undulationFilename;

% Requested to use "high resolution", use hi-res geoid file regardless of
% the passed-in filename parameter.
if p.Results.useHiRes
    filename = defaultHiResUundulationFilename;
end


% Note: we need to append the '.' to the file name because "which" uses
% that to identify the string as a file (rather than a function).  If
% the file name already has a '.' it still works, but if it's missing
% one "which" will fail.
fullPath = which([filename '.']);

fileStats = dir(fullPath);
if isempty(fileStats),
    error('GEOID_UNDULATION:UNDULATION_FILE_NOT_FOUND',...
        'Undulation file not found.  Pass filename explicitly with the ''undulationFilename'' input parameter, or place a file with the default filename on the MATLAB path.');
end;
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


fid = fopen(filename, 'r');

% Seek to the correct location in the file.  Each record is a set
% of longitudes at a constant latitude.  Seek to the record which is
% just "above" the desired latitude.  If -90 degrees is being requested
% grab the last two records.  We'll do that by "cheating" and forcing
% the largest offset value to be one less than the number of records.  
% later we'll grab two records so we'll get the last one.
off = min( floor((90-lat)/delta_lat), num_records-1 );
fseek(fid, off*record_length, 'bof');

% Read in two records.  This gives us the latitude "cuts" above and
% below the desired latitude value.  We'll also strip out the leading
% and trailing record markers and extend the values to "wrap".  The
% file doesn't duplicate the value at 360 degrees (which should equal
% the value at 0 degreees) so we'll add it in - that is, copy the value
% from 0 degrees to the end of the lon_vals array (360 degrees).
lon_vals = fread(fid, [num_record_entries 2], 'float32')';
fclose(fid);
lon_vals = lon_vals(:,2:end-1);
lon_vals(:,num_lon_values+1) = lon_vals(:,1);

% Now pull out the point values surrounding the requested longitude.
% Here "index" is the index of the point nearest west so the
% requested longitude must lie between index and index+1.
index = floor(1 + lon/delta_lon);
vals = lon_vals(:,index:index+1);

% We'll use Matlab's built-in interpolator so we'll need to build
% arrays to represet the latitude and longitude "grid" of the four
% points.  The "Xs" and "Ys" arrays give the latitudes and longitudes 
% of the four corners while "vals" gives the undulation at those
% points.
Xs = delta_lon*[index-1 index;index-1 index];
Ys = 90 - delta_lat*[off off;off+1 off+1];    
undulation = interp2(Xs, Ys, vals, lon, lat, 'linear');
    
end
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

