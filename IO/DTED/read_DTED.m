function [elevations, lats, lons, meta] = read_DTED( filename, LL, UR )
%read_DTED Read DTED formatted elevation data
%   [elevations, lats, lons, meta] = read_DTED( filename )
%
% This reads files in the format described in:
%      Performance Specification Digital Terrain Elevation Data (DTED)
%      MIL-PRF-89020B
%      23 May 2000
%
%   INPUTS:
%      FILENAME:      Input file to read (DTED format).
%      LL:            A two-vector of latitude, longitude (in degrees) of 
%                     the lower-left corner of the region to load.
%      UR:            A two-vector of latitude, longitude (in degrees) of
%                     the upper-right corner of the region to load.
%
%   OUTPUTS:
%      ELEVATIONS:    A matrix of elevation values in meters above 
%                     mean sea level.  Each row is a "slice" of 
%                     constant-latitude values in increasing latitude.
%                     That is, point (1,1) is the smallest latitude,
%                     smallest longitude.  Point (N,N) is highest latitude
%                     highest longitude.
%      LATS:          The latitude values corresponding to each row of
%                     ELEVATIONS.  These are in decimal degrees (double).
%      LONS:          The longitude values corresponding to each row of
%                     ELEVATIONS.  These are in decimal degrees (double).
%      META:          Supporting meta data.  This is a structure holding
%                     _some_ (not all) of the data from the various header
%                     records (user, data set, accuracy, etc.)
%
%   Author: Tom Krauss (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

    fid = fopen(filename);
    if fid == -1
        error('read_DTED:nofile', 'Unable to open file "%s"\n', filename);
    end

    meta = struct();
    
    % Check to make sure we have a user header for an uncompressed DTED
    % file (UHL in first three bytes).
    uhlhdr = fread(fid, 3, 'uint8=>char')';
    if strcmp(uhlhdr, 'UHL') ~= 1
        error('read_DTED:notDTED', 'File not DTED format or unsupported (compressed) DTED format');
    end

    % Jump over the rest of the user header and go right to the 
    % Data Set Identifier record.
    fseek(fid, 80, 'bof');
    dsihdr    = fread(fid, 3, 'uint8=>char')';
    if strcmp(dsihdr, 'DSI') ~= 1
        error('read_DTED:badFile', 'Data Set Identifier Record not found: File not DTED format?');
    end
    meta.sec_class = fread(fid, 1, 'uint8=>char')';
    meta.sec_mark  = fread(fid, 2, 'uint8=>char')';
    meta.sec_hand  = fread(fid, 27, 'uint8=>char')';
    fread(fid, 26, 'uint8=>char'); % fill

    meta.dma_desc = fread(fid, 5, 'uint8=>char')';
    fread(fid, 23, 'uint8=>char'); % fill

    meta.data_edition    = fread(fid, 2, 'uint8=>char')';
    meta.match_version   = fread(fid, 1, 'uint8=>char')';
    meta.maint_date      = fread(fid, 4, 'uint8=>char')';
    meta.match_date      = fread(fid, 4, 'uint8=>char')';
    meta.maint_desc_code = fread(fid, 4, 'uint8=>char')';
    meta.producer_code   = fread(fid, 8, 'uint8=>char')';
    fread(fid, 16, 'uint8=>char'); % fill

    meta.product_vers = fread(fid, 9, 'uint8=>char')';
    fread(fid, 2, 'uint8=>char'); % fill

    meta.product_date = fread(fid, 4, 'uint8=>char')';
    meta.vert_datum   = fread(fid, 3, 'uint8=>char')';
    meta.horiz_datum  = fread(fid, 5, 'uint8=>char')';
    meta.dig_coll_sys = fread(fid, 10, 'uint8=>char')';
    meta.comp_date    = fread(fid, 4, 'uint8=>char')';
    fread(fid, 22, 'uint8=>char'); % fill

    lat_orig_deg     = str2double(fread(fid, 2, 'uint8=>char')');
    lat_orig_min     = str2double(fread(fid, 2, 'uint8=>char')');
    lat_orig_sec     = str2double(fread(fid, 4, 'uint8=>char')');
    lat_orig_dir     = fread(fid, 1, 'uint8=>char')';
    if (lat_orig_dir == 'S'  ||  lat_orig_dir == 's'  || ...
        lat_orig_dir == 'W'  ||  lat_orig_dir == 'w')
        lat_orig_deg = -lat_orig_deg;
    end
    meta.lat_orig_dd = sign(lat_orig_deg)*polyval(abs([lat_orig_sec lat_orig_min lat_orig_deg]),1/60);
    
    lon_orig_deg     = str2double(fread(fid, 3, 'uint8=>char')');
    lon_orig_min     = str2double(fread(fid, 3, 'uint8=>char')');
    lon_orig_sec     = str2double(fread(fid, 3, 'uint8=>char')');
    lon_orig_dir     = fread(fid, 1, 'uint8=>char')';
    if (lon_orig_dir == 'S'  ||  lon_orig_dir == 's'  || ...
        lon_orig_dir == 'W'  ||  lon_orig_dir == 'w')
        lon_orig_deg = -lon_orig_deg;
    end
    meta.lon_orig_dd = sign(lat_orig_deg)*polyval(abs([lon_orig_sec lon_orig_min lon_orig_deg]),1/60);
    
    meta.lat_sw_corn      = fread(fid, 7, 'uint8=>char')';
    meta.lon_sw_corn      = fread(fid, 8, 'uint8=>char')';
    meta.lat_nw_corn      = fread(fid, 7, 'uint8=>char')';
    meta.lon_nw_corn      = fread(fid, 8, 'uint8=>char')';
    meta.lat_ne_corn      = fread(fid, 7, 'uint8=>char')';
    meta.lon_ne_corn      = fread(fid, 8, 'uint8=>char')';
    meta.lat_se_corn      = fread(fid, 7, 'uint8=>char')';
    meta.lon_se_corn      = fread(fid, 8, 'uint8=>char')';
    meta.orient_angle     = str2double(fread(fid, 9, 'uint8=>char')');
    meta.lat_spacing_dd   = str2double(fread(fid, 4, 'uint8=>char')')/(10*60*60);
    meta.lon_spacing_dd   = str2double(fread(fid, 4, 'uint8=>char')')/(10*60*60);
    meta.num_lat_lines    = str2double(fread(fid, 4, 'uint8=>char')');
    meta.num_lon_lines    = str2double(fread(fid, 4, 'uint8=>char')');
    meta.part_cell_ind    = fread(fid, 2, 'uint8=>char')';
    fread(fid, 357, 'uint8=>char'); % fill

    % Accuracy record...
    accdes         = fread(fid, 3, 'uint8=>char')';
    if strcmp(accdes, 'ACC') ~= 1
        error('read_DTED:badFile', 'Accuracy Description Record not found: File not DTED format?');
    end
    meta.abs_hor_acc    = str2double(fread(fid, 4, 'uint8=>char')');
    meta.abs_vert_acc   = str2double(fread(fid, 4, 'uint8=>char')');
    meta.pt2pt_hor_acc  = str2double(fread(fid, 4, 'uint8=>char')');
    meta.pt2pt_vert_acc = str2double(fread(fid, 4, 'uint8=>char')');
    fread(fid, 36, 'uint8=>char'); % fill

    meta.mult_acc_out_flg = fread(fid, 2, 'uint8=>char')';

    % Build the reference latitude and longitude values
    r = meta.num_lat_lines;
    c = meta.num_lon_lines;
    %PJC 2014-05-21 fixed bug for negative lats or lons
    lats = (0:(c-1))*meta.lat_spacing_dd + meta.lat_orig_dd*sign(LL(1));
    lons = (0:(r-1))*meta.lon_spacing_dd + meta.lon_orig_dd*sign(LL(2));


    % If range of lat/lons is given, constrain window to read
    if nargin>2
        latTmp = [find(lats<=LL(1),1,'last') find(lats<=UR(1),1,'last')];
        lower_lat_index = min(latTmp);
        upper_lat_index = max(latTmp);
        lats = lats(lower_lat_index:upper_lat_index);
        
        lonTmp = [find(lons<=LL(2),1,'last') find(lons<=UR(2),1,'last')];
        lower_lon_index = min(lonTmp);
        upper_lon_index = max(lonTmp);
        lons = lons(lower_lon_index:upper_lon_index);
    else
        lower_lat_index=1;
        upper_lat_index=meta.num_lat_lines;
        lower_lon_index=1;
        upper_lon_index=meta.num_lon_lines;
    end
    
    % Skip to the correct record and read data (the bracketing longitudes),
    % then extract the bracketing latitudes.
    %
    % Each row of bulk_data is a record.  Each record is composed of:
    %      1 byte sentinal (decimal value 170)
    %      3 byte record count
    %      2 byte longitude count
    %      2 byte latitude count
    %      N*2 byte elevation values
    %      4 byte checksum
    data_record_length = 1 + 3 + 2 + 2 + 2*meta.num_lat_lines + 4;
    fseek(fid, 3429-1 + (lower_lon_index-1)*data_record_length, 'bof');
    bulk_data = fread(fid, [data_record_length (upper_lon_index-lower_lon_index+1)], 'uint8')';
    fclose(fid);

    % Check all record sentinels (first column)
    check = find( (bulk_data(:,1) ~= 170)   &  (bulk_data(:,1) ~= 0), 1 );
    if ~isempty(check)
        error('read_DTED:badFile', 'Record sentinal not found.  File not in DTED format?');
    end

    % Extract out the "high" and "low" bytes
    high = bulk_data(:,9:2:data_record_length-4);
    low  = bulk_data(:,10:2:data_record_length-4);
    if (~ all (size(high)==size(low)) )
        error('read_DTED:badFile', 'ERROR: Data records wrong size(?)');
    end

    % Check the sign bit on the "high" byte.  If it's set, unset it and change
    % the sign of the byte.
    temp_elev = ones(size(high));
    idx = find( high>128 );
    high(idx) = high(idx)-128;
    temp_elev(idx) = -1;

    % Make the actual elevation word (shift high and add)
    elevations = (temp_elev(:,lower_lat_index:upper_lat_index) .* ...
        (high(:,lower_lat_index:upper_lat_index)*256 + ...
        low(:,lower_lat_index:upper_lat_index)));

    % Clean up voids
    idx = elevations < -50;
    elevations(idx) = -50;
end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////