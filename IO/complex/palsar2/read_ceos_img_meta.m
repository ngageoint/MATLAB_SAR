function [ meta ] = read_ceos_img_meta( filename )
%READ_CEOS_IMG_META Read CEOS SAR image file
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open file
fid=fopen(filename,'r','b');

% Read volume descriptor records
meta.rec_num=fread(fid,1,'uint32'); % Record sequence number
meta.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
meta.rec_type=fread(fid,1,'uint8'); % Record type
meta.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
meta.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
meta.rec_length=fread(fid,1,'uint32'); % Record length
meta.ascii_ebcdic=fread(fid,2,'uint8=>char').'; % Record length
fseek(fid,2,'cof'); % Blanks
meta.doc_id=fread(fid,12,'uint8=>char').'; % Superstructure format control document ID
meta.doc_rev=fread(fid,2,'uint8=>char').'; % Superstructure format control document revision level
meta.rec_rev=fread(fid,2,'uint8=>char').'; % Superstructure record format revision level
meta.soft_rel_rev=fread(fid,12,'uint8=>char').'; % Software release and revision level
meta.file_num=fread(fid,4,'uint8=>char').'; % File number
meta.file_id=fread(fid,16,'uint8=>char').'; % File ID
meta.rec_seq_loc_type_flag=fread(fid,4,'uint8=>char').'; % Record sequency and location type flag
meta.seq_num_loc=fread(fid,8,'uint8=>char').'; % Sequence number of location
meta.fld_len_seq=fread(fid,4,'uint8=>char').'; % Field length of sequence number
meta.rec_code_loc_type_flag=fread(fid,4,'uint8=>char').'; % Records code and location type flag
meta.loc_rec_code=fread(fid,8,'uint8=>char').'; % Location of record code
meta.fld_len_code=fread(fid,4,'uint8=>char').'; % Field length of record code
meta.rec_len_loc_type_flag=fread(fid,4,'uint8=>char').'; % Record length and location type flag
meta.loc_rec_len=fread(fid,8,'uint8=>char').'; % Location of record length
meta.len_rec_len=fread(fid,4,'uint8=>char').'; % Location of record length
fseek(fid,68,'cof'); % Blanks
meta.num_data_rec=str2double(fread(fid,6,'uint8=>char')); % Number of data set records
meta.data_len=str2double(fread(fid,6,'uint8=>char')); % Data set summary record length
fseek(fid,24,'cof'); % Blanks
% Sample group data
meta.sample_len=str2double(fread(fid,4,'uint8=>char')); % Bit length per sample
meta.num_samples=str2double(fread(fid,4,'uint8=>char')); % Number of samples per data group
meta.num_bytes=str2double(fread(fid,4,'uint8=>char')); % Number of bytes per data group
meta.just_order=fread(fid,4,'uint8=>char').'; % Justification and order of samples
% SAR related data
meta.num_chan=str2double(fread(fid,4,'uint8=>char')); % Number of SAR channels
meta.num_lines=str2double(fread(fid,8,'uint8=>char')); % Number of lines per data set
meta.num_left=str2double(fread(fid,4,'uint8=>char')); % Number of left border pixels per line
meta.num_pixels=str2double(fread(fid,8,'uint8=>char')); % Number of data groups (or pixels) per line
meta.num_right=str2double(fread(fid,4,'uint8=>char')); % Number of right border pixels per line
meta.num_top=str2double(fread(fid,4,'uint8=>char')); % Number of top border lines
meta.num_bottom=str2double(fread(fid,4,'uint8=>char')); % Number of bottom border lines
meta.interleave=fread(fid,4,'uint8=>char').'; % Interleaving ID
% Record data
meta.phys_rec_line=str2double(fread(fid,2,'uint8=>char')); % Number of physical records per line
meta.phys_rec_multi_chan=str2double(fread(fid,2,'uint8=>char')); % Number of physical records per multi-channel line
meta.prefix_bytes=str2double(fread(fid,4,'uint8=>char')); % Number of bytes of prefix per record
meta.sar_data_bytes=str2double(fread(fid,8,'uint8=>char')); % Number of bytes of SAR data per record
meta.suffix_bytes=str2double(fread(fid,4,'uint8=>char')); % Number of bytes of suffix data per record
meta.pre_suf_rpt_flg=fread(fid,4,'uint8=>char').'; % Prefix/suffix repeat flag
% Prefix/suffix data locations
meta.loc_sar_data=fread(fid,8,'uint8=>char').'; % Sample data line number locator
meta.loc_sar_chan_num=fread(fid,8,'uint8=>char').'; % SAR channel number locator
meta.loc_time=fread(fid,8,'uint8=>char').'; % Time of SAR data line locator
meta.loc_leftfill=fread(fid,8,'uint8=>char').'; % Left-fill count locator
meta.loc_rightfill=fread(fid,8,'uint8=>char').'; % Right-fill count locator
meta.pad_pixels=fread(fid,4,'uint8=>char').'; % Pad pixels present indictor
fseek(fid,28,'cof'); % Blanks
meta.loc_data_qual=fread(fid,8,'uint8=>char').'; % SAR data line quality code locator
meta.loc_cal_info=fread(fid,8,'uint8=>char').'; % Calibration information field locator
meta.loc_gain=fread(fid,8,'uint8=>char').'; % Gain values field locator
meta.loc_bias=fread(fid,8,'uint8=>char').'; % Bias values field lcoator
meta.sar_datatype=fread(fid,28,'uint8=>char').'; % SAR data format type indicator
meta.sar_datatype_code=fread(fid,4,'uint8=>char').'; % SAR data format type code
meta.num_leftfill=fread(fid,4,'uint8=>char').'; % Number of left fill bits within pixel
meta.num_rightfill=fread(fid,4,'uint8=>char').'; % Number of right fill bits within pixel
meta.max_data_range=fread(fid,8,'uint8=>char').'; % Maximum data range of pixel
meta.scansar_num_bursts=fread(fid,4,'uint8=>char').'; % ScanSAR, number of burst data in this file
meta.scansar_num_lines=fread(fid,4,'uint8=>char').'; % ScanSAR, number of lines per one burst
meta.scansar_num_overlap=fread(fid,4,'uint8=>char').'; % ScanSAR, number of overlap lines with adjacent bursts
fseek(fid,260,'cof'); % Blanks
sig_start = ftell(fid);

% Pixel data
switch meta.file_id(8)
    case 'B'
        % Signal data records
        % recs = 1:meta.num_data_rec;
        % Only need the first (and maybe last) record, since these
        % per-record fields are mostly all identical except for a few
        % fields we don't need.
        recs = [1, meta.num_data_rec];
        for i = 1:numel(recs)
            fseek(fid, sig_start + (meta.prefix_bytes + ...
                (meta.num_pixels*meta.num_bytes)) * (recs(i)-1), 'bof');
            meta.signal(i).rec_num=fread(fid,1,'uint32'); % Record sequence number
            meta.signal(i).rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
            meta.signal(i).rec_type=fread(fid,1,'uint8'); % Record type
            meta.signal(i).rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
            meta.signal(i).rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
            meta.signal(i).rec_length=fread(fid,1,'uint32'); % Record length
            % Prefix data-general information
            meta.signal(i).line_num=fread(fid,1,'uint32'); % SAR image data line number
            meta.signal(i).sar_rec_ind=fread(fid,1,'uint32'); % SAR image data record index
            meta.signal(i).left_fill=fread(fid,1,'uint32'); % Actual count of left-fill pixels
            meta.signal(i).num_pixels=fread(fid,1,'uint32'); % Actual count of data pixels
            meta.signal(i).right_fill=fread(fid,1,'uint32'); % Actual count of right-fill pixels
            % Prefix data-sensor parameters
            meta.signal(i).update_flg=fread(fid,1,'uint32'); % Sensor parameters update flag
            meta.signal(i).year=fread(fid,1,'uint32'); % Sensor acquisition year
            meta.signal(i).day=fread(fid,1,'uint32'); % Sensor acquisition day of year
            meta.signal(i).msec=fread(fid,1,'uint32'); % Sensor acquisition milli-seconds of day
            meta.signal(i).chan_id=fread(fid,1,'uint16'); % SAR channel ID
            meta.signal(i).chan_code=fread(fid,1,'uint16'); % SAR channel code
            meta.signal(i).tx_pol=fread(fid,1,'uint16'); % Transmitted polarization
            meta.signal(i).rcv_pol=fread(fid,1,'uint16'); % Received polarization
            meta.signal(i).prf=fread(fid,1,'uint32'); % PRF [mHz]
            meta.signal(i).scan_id=fread(fid,1,'uint32'); % Scan ID
            meta.signal(i).rng_comp_flg=fread(fid,1,'uint16'); % Onboard range compressed flag
            meta.signal(i).chirp_type=fread(fid,1,'uint16'); % Chirp type designator
            meta.signal(i).chirp_length=fread(fid,1,'uint32'); % Chirp length (pulse width) [nsec]
            meta.signal(i).chirp_const=fread(fid,1,'int32'); % Chirp constant coefficient [Hz]
            meta.signal(i).chirp_lin=fread(fid,1,'int32'); % Chirp linear coefficient [Hz/usec]
            meta.signal(i).chirp_quad=fread(fid,1,'int32'); % Chirp quadratic coefficient [Hz/usec^2]
            meta.signal(i).usec=fread(fid,1,'uint64'); % Sensor acquisition micro-second of day (rounded or floored)
            meta.signal(i).gain=fread(fid,1,'uint32'); % Receiver gain [dB]
            meta.signal(i).invalid_flg=fread(fid,1,'uint32'); % Invalid line flag
            meta.signal(i).elec_ele=fread(fid,1,'int32'); % Electronic elevation angle at nadir of antenna [deg]
            meta.signal(i).mech_ele=fread(fid,1,'int32'); % Mechanical elevation angle at nadir of antenna [deg]
            meta.signal(i).elec_squint=fread(fid,1,'int32'); % Electronic antenna squint angle [deg]
            meta.signal(i).mech_squint=fread(fid,1,'int32'); % Mechanical antenna squint angle [deg]
            meta.signal(i).slant_rng=fread(fid,1,'uint32'); % Slant range to 1st data sample [m]
            meta.signal(i).wind_pos=fread(fid,1,'uint32'); % Data record window position (SAMPLE DELAY [nsec])
            fseek(fid,4,'cof'); % Blanks
            % Prefix data-platform reference information
            meta.signal(i).pos_update_flg=fread(fid,1,'uint32'); % Platform position parameters update flag
            meta.signal(i).plat_lat=fread(fid,1,'int32'); % Platform latitude [1/1,000,000 deg]
            meta.signal(i).plat_lon=fread(fid,1,'int32'); % Platform longitude [1/1,000,000 deg]
            meta.signal(i).plat_alt=fread(fid,1,'int32'); % Platform altitude [m]
            meta.signal(i).grnd_spd=fread(fid,1,'int32'); % Platform ground speed [cm/sec]
            meta.signal(i).vel_x=fread(fid,1,'int32'); % Platform velocity X [cm/sec]
            meta.signal(i).vel_y=fread(fid,1,'int32'); % Platform velocity Y [cm/sec]
            meta.signal(i).vel_z=fread(fid,1,'int32'); % Platform velocity Z [cm/sec]
            meta.signal(i).acc_x=fread(fid,1,'int32'); % Platform acceleration X [cm/sec^2]
            meta.signal(i).acc_y=fread(fid,1,'int32'); % Platform acceleration X [cm/sec^2]
            meta.signal(i).acc_z=fread(fid,1,'int32'); % Platform acceleration X [cm/sec^2]
            meta.signal(i).track=fread(fid,1,'int32'); % Platform track angle [1/1,000,000 deg]
            meta.signal(i).true_track=fread(fid,1,'int32'); % Platform rue track angle [1/1,000,000 deg]
            meta.signal(i).pitch=fread(fid,1,'int32'); % Platform pitch angle [1/1,000,000 deg]
            meta.signal(i).roll=fread(fid,1,'int32'); % Platform roll angle [1/1,000,000 deg]
            meta.signal(i).yaw=fread(fid,1,'int32'); % Platform yaw angle [1/1,000,000 deg]
            % Prefix data-sensor/facility specific auxiliary data
            meta.signal(i).lat_first=fread(fid,1,'int32'); % Latitude of 1st pixel [1/1,000,000 deg]
            meta.signal(i).lat_center=fread(fid,1,'int32'); % Latitude of center pixel [1/1,000,000 deg]
            meta.signal(i).lat_last=fread(fid,1,'int32'); % Latitude of last pixel [1/1,000,000 deg]
            meta.signal(i).lon_first=fread(fid,1,'int32'); % Longitude of 1st pixel [1/1,000,000 deg]
            meta.signal(i).lon_center=fread(fid,1,'int32'); % Longitude of center pixel [1/1,000,000 deg]
            meta.signal(i).lon_last=fread(fid,1,'int32'); % Longitude of last pixel [1/1,000,000 deg]
            % ScanSAR burst data parameters
            meta.signal(i).burst_num=fread(fid,1,'uint32'); % ScanSAR, burst number
            meta.signal(i).line_num=fread(fid,1,'uint32'); % ScanSAR, line number in this burst
            fseek(fid,60,'cof'); % Blanks
            meta.signal(i).frame_num=fread(fid,1,'uint32'); % ALOS2 frame number
            % fseek(fid,256,'cof'); % PALSAR auxiliary data
            % fseek(fid,meta.num_pixels*meta.num_bytes,'cof'); % SAR data
        end
    case {'C','D'}
        % Processed data records
        disp('Processed data');
end

% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////