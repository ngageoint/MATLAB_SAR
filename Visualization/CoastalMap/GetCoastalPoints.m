function [lats, lons] = GetCoastalPoints(Quality,maxLevel,minArea,LatLon)
%GETCOASTALPOINTS Reads the coastal points from a GSHHS file (version 2.x)
%
% USAGE:
%   [lats, lons] = GetCoastalPoints(Quality,maxLevel,minArea,LatLon)
%
% INPUTS:   
%   Quality - integer specifying the quality of the coastline to load.
%       Acceptable values are 1-5. Only the coastline corresponding to 1
%       is available by default. REQUIRED
%   maxLevel - max level of data coastline to load. Levels are:
%           1 land, 2 lake, 3 island in lake, 4 pond in island in lake.
%           Default 1
%   minArea - minimum area (km^2) within a coastline to load. Default 10
%   LatLon - 1 by 4 vector describing a window in which a coastline must be
%       included in order for it to be loaded. Values are 
%       [minLat minLong maxLat maxLon]. If empty load all. Default []
%
% OUTPUTS:
%   lats - latitude coordinates.
%   lons - longitude coordinates.
%
% Note: shapes are separated by NaN
%
% GSHHS data and a description of its format can be found at:
%
% "A Global Self-consistent, Hierarchical, High-resolution Shoreline Database"
% http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
%
% or
%
% "A Global Self-consistent, Hierarchical, High-resolution Geography Database"
% http://www.soest.hawaii.edu/pwessel/gshhs/index.html
%
% Author: Sean Hatch
% Modified by Wade Schwartzkopf to handle GSHHS version 2.x
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%Modified by Pat Cutler (NGA) 20130306 for load images segments LatLon
%input and documented inputs and outputs

if ~exist('maxLevel','var') || isempty(maxLevel)
    maxLevel = 1;
end

if ~exist('minArea','var') || isempty(minArea)
    minArea = 10;
end

if ~exist('LatLon','var') || isempty(LatLon)
    LatLonCheck = 0;
else
    LatLonCheck = 1;
end

% open the GSHHS data file
% the _c file contains the crudest data in the shoreline hierarchy
files = {'gshhs_c.b', 'gshhs_l.b', 'gshhs_i.b', 'gshhs_h.b', 'gshhs_f.b'};
coastline_file = files{Quality};
fid = fopen(coastline_file, 'r', 'b');
% read in all of the data blocks from the file
tot_points = 1;
lons = NaN;
lats = NaN;
while true
    % read the header of the block
    [polygon_id, count] = fread(fid, 1, 'int32'); % Unique polygon id number, starting at 0
    if (count < 1)
        break;
    end
    num_points = fread(fid, 1, 'int32');  % Number of points in this polygon
    flag = fread(fid, 1, 'int32');
    level = bitand(flag, 255); % 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
    version = bitand(bitshift(flag,-8),255); % Should be 7 for GSHHS release 7 (i.e., version 2.0)
    greenwich = bitget(flag,17); % Greenwich is 1 if Greenwich is crossed
    source = bitget(flag,25); % 0 = CIA WDBII, 1 = WVS
    river = bitget(flag,26); % 0 = not set, 1 = river-lake and level = 2
    % Min/max extent stored in file in micro-degrees
    west  = fread(fid, 1, 'int32') ./ 1e6; % degrees
    east  = fread(fid, 1, 'int32') ./ 1e6; % degrees
    south = fread(fid, 1, 'int32') ./ 1e6; % degrees
    north = fread(fid, 1, 'int32') ./ 1e6; % degrees
    area  = fread(fid, 1, 'int32') ./10;   % Area of polygon in km^2
    area_full  = fread(fid, 1, 'int32') ./10; % Area of original full-resolution polygon in km^2
    container = fread(fid, 1, 'int32'); % Id of container polygon that encloses this polygon (-1 if none)
    % Id of ancestor polygon in the full resolution set that was the source of this polygon (-1 if none)
    ancestor = fread(fid, 1, 'int32');
    % make sure that this is a coast for land of a reasonable size
    if (level <= maxLevel) && (area(1) > minArea)
        if checkAOI(east,west,north,south)
            % read in the lons/lats of all the data points
            data = fread(fid, [2,num_points] , 'int32') ./ 1e6;
            
            % shift the range of the longitudes
            hi = find(data(1,:) > 180);
            data(1,hi) = data(1,hi) - 360;
            
            % check for wrap around
            wraps = [find(abs(diff(data(1, :))) > 180), num_points];
            num_wraps = length(wraps);
            unwrap = data(:,1:wraps(1));
            if num_wraps > 1
                check = zeros(num_wraps-1,1);
                for index = 2:num_wraps
                    check(index-1) = checkAOI(max(data(2,wraps(index-1)+1:wraps(index))),...
                        min(data(2,wraps(index-1)+1:wraps(index))),...
                        max(data(1,wraps(index-1)+1:wraps(index))),...
                        min(data(1,wraps(index-1)+1:wraps(index))));
                    
                    unwrap = [unwrap, [NaN;NaN], data(:,wraps(index-1)+1:wraps(index))];
                end
                if ~any(check)
                    continue
                end
            end
            data = unwrap;
            num_points = num_points + num_wraps - 1;
            
            % concatenate them to the rest of the list
            tot_points = tot_points + num_points + 1;
            lons = [lons, data(1, :), NaN];
            lats = [lats, data(2, :), NaN];
        else
            % skip over the data points
            fseek(fid, 8.*num_points, 'cof');
        end
    else
        % skip over the data points
        fseek(fid, 8.*num_points, 'cof');
    end
end

% close the data file
fclose(fid);


    function check = checkAOI(east,west,north,south)
        if LatLonCheck
            if west > 180
                west = west-360;
            end
            if east > 180
                east = east-360;
            end
            %longitudinal exclusion
            LonExclude = (west < LatLon(2) && east < LatLon(2)) || (west > LatLon(4) && east > LatLon(4));
            %latitudinal exclusion
            LatExclude = (south < LatLon(1) && north < LatLon(1)) || (south > LatLon(3) && north > LatLon(3));
            if LonExclude || LatExclude
                check = 0;
            else
                %longitudinal crossing
                LonCross = (west < LatLon(2) && east > LatLon(2)) || (west < LatLon(4) && east > LatLon(4));
                %latitudinal crossing
                LatCross = (south < LatLon(1) && north > LatLon(1)) || (south < LatLon(3) && north > LatLon(3));
                if LonCross && ~LonExclude
                    check = 1;
                else if LatCross && ~LatExclude
                        check = 1;
                    else
                        check = 0;
                    end
                end
            end
        else
            check = 1;
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////