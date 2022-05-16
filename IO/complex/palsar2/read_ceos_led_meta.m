function [ meta ] = read_ceos_led_meta( filename )
%READ_CEOS_LED_META Read CEOS SAR leader file
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open file
fid=fopen(filename,'r','b');

% Read volume descriptor records
meta.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
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
meta.num_map_rec=str2double(fread(fid,6,'uint8=>char')); % Number of map projection data records
meta.map_len=str2double(fread(fid,6,'uint8=>char')); % Map projection data record length
meta.num_pos_rec=str2double(fread(fid,6,'uint8=>char')); % Number of platform position data records
meta.pos_len=str2double(fread(fid,6,'uint8=>char')); % Platform position record length
meta.num_att_rec=str2double(fread(fid,6,'uint8=>char')); % Number of attitude data records
meta.att_len=str2double(fread(fid,6,'uint8=>char')); % Attitude data record length
meta.num_rad_rec=str2double(fread(fid,6,'uint8=>char')); % Number of radiometric data records
meta.rad_len=str2double(fread(fid,6,'uint8=>char')); % Radiometric data record length
meta.num_rad_comp_rec=str2double(fread(fid,6,'uint8=>char')); % Number of radiometric compensation records
meta.rad_comp_len=str2double(fread(fid,6,'uint8=>char')); % Radiometric compensation record length
meta.num_data_qual_rec=str2double(fread(fid,6,'uint8=>char')); % Number of data quality summary records
meta.data_qual_len=str2double(fread(fid,6,'uint8=>char')); % Data quality summary record length
meta.num_hist_rec=str2double(fread(fid,6,'uint8=>char')); % Number of data histogram records
meta.hist_len=str2double(fread(fid,6,'uint8=>char')); % Data histogram record length
meta.num_rng_spect_rec=str2double(fread(fid,6,'uint8=>char')); % Number of range spectra records
meta.rng_spect_len=str2double(fread(fid,6,'uint8=>char')); % Range spectra record length
meta.num_dem_rec=str2double(fread(fid,6,'uint8=>char')); % Number of DEM descriptor records
meta.dem_len=str2double(fread(fid,6,'uint8=>char')); % DEM descriptor record length
meta.num_radar_rec=str2double(fread(fid,6,'uint8=>char')); % Number of radar parameter update records
meta.radar_len=str2double(fread(fid,6,'uint8=>char')); % Radar parameter update record length
meta.num_annot_rec=str2double(fread(fid,6,'uint8=>char')); % Number of annotation data records
meta.annot_len=str2double(fread(fid,6,'uint8=>char')); % Annotation data record length
meta.num_proc_rec=str2double(fread(fid,6,'uint8=>char')); % Number of detail processing records
meta.proc_len=str2double(fread(fid,6,'uint8=>char')); % Detail processing record length
meta.num_cal_rec=str2double(fread(fid,6,'uint8=>char')); % Number of calibration records
meta.cal_len=str2double(fread(fid,6,'uint8=>char')); % Calibration record length
meta.num_gcp_rec=str2double(fread(fid,6,'uint8=>char')); % Number of GCP records
meta.gcp_len=str2double(fread(fid,6,'uint8=>char')); % GCP record length
fseek(fid,60,'cof'); % Blanks
meta.num_fac_data_rec1=str2double(fread(fid,6,'uint8=>char')); % Number of facility data records
meta.fac_data_len1=str2double(fread(fid,8,'uint8=>char')); % Facility data record length
meta.num_fac_data_rec2=str2double(fread(fid,6,'uint8=>char')); % Number of facility data records
meta.fac_data_len2=str2double(fread(fid,8,'uint8=>char')); % Facility data record length
meta.num_fac_data_rec3=str2double(fread(fid,6,'uint8=>char')); % Number of facility data records
meta.fac_data_len3=str2double(fread(fid,8,'uint8=>char')); % Facility data record length
meta.num_fac_data_rec4=str2double(fread(fid,6,'uint8=>char')); % Number of facility data records
meta.fac_data_len4=str2double(fread(fid,8,'uint8=>char')); % Facility data record length
meta.num_fac_data_rec5=str2double(fread(fid,6,'uint8=>char')); % Number of facility data records
meta.fac_data_len5=str2double(fread(fid,8,'uint8=>char')); % Facility data record length
fseek(fid,230,'cof'); % Blanks

% Read data set summary
meta.data.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
meta.data.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
meta.data.rec_type=fread(fid,1,'uint8'); % Record type
meta.data.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
meta.data.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
meta.data.rec_length=fread(fid,1,'uint32'); % Record length
meta.data.rec_seq=fread(fid,4,'uint8=>char').'; % Record sequence number
meta.data.sar_id=fread(fid,4,'uint8=>char').'; % SAR channel ID
meta.data.scene_id=fread(fid,32,'uint8=>char').'; % Scene ID
meta.data.num_scene_ref=fread(fid,16,'uint8=>char').'; % Number of scene reference
meta.data.scene_ctr_time=fread(fid,32,'uint8=>char').'; % Scene center time
fseek(fid,16,'cof'); % Blanks
meta.data.geo_lat=fread(fid,16,'uint8=>char').'; % Geodetic latitude of processed scene center
meta.data.geo_long=fread(fid,16,'uint8=>char').'; % Geodetic longitude of processed scene center
meta.data.heading=fread(fid,16,'uint8=>char').'; % Processed scene center true heading
meta.data.ellips=fread(fid,16,'uint8=>char').'; % Ellipsoid designator
meta.data.semimajor=fread(fid,16,'uint8=>char').'; % Ellipsoid semi-major axis
meta.data.semiminor=fread(fid,16,'uint8=>char').'; % Ellipsoid semi-minor axis
meta.data.earth_mass=fread(fid,16,'uint8=>char').'; % Earth mass
meta.data.grav_const=fread(fid,16,'uint8=>char').'; % Gravitational constant
meta.data.J2=fread(fid,16,'uint8=>char').'; % Ellipsoid J2 parameter
meta.data.J3=fread(fid,16,'uint8=>char').'; % Ellipsoid J3 parameter
meta.data.J4=fread(fid,16,'uint8=>char').'; % Ellipsoid J4 parameter
fseek(fid,16,'cof'); % Blanks
meta.data.avg_terr=fread(fid,16,'uint8=>char').'; % Average terrain height abovie ellipsoid at scene center
meta.data.ctr_line=str2double(fread(fid,8,'uint8=>char')); % Scene center line number
meta.data.ctr_pixel=str2double(fread(fid,8,'uint8=>char')); % Scene center pixel number
meta.data.proc_length=fread(fid,16,'uint8=>char').'; % Processing scene length [km]
meta.data.proc_width=fread(fid,16,'uint8=>char').'; % Processing scene width [km]
fseek(fid,16,'cof'); % Blanks
meta.data.sar_chan=fread(fid,4,'uint8=>char').'; % Number of SAR channel
fseek(fid,4,'cof'); % Blanks
meta.data.platform=fread(fid,16,'uint8=>char').'; % Sensor platform mission identifier
meta.data.sensor_id_mode=fread(fid,32,'uint8=>char').'; % Sensor ID and operational mode
meta.data.orbit=fread(fid,8,'uint8=>char').'; % Orbit number or flight line
meta.data.sensor_lat=fread(fid,8,'uint8=>char').'; % Sensor platform geodetic latitude at nadir corresponding to scene center
meta.data.sensor_lon=fread(fid,8,'uint8=>char').'; % Sensor platform geodetic longitude at nadir corresponding to scene center
meta.data.sensor_heading=fread(fid,8,'uint8=>char').'; % Sensor platform heading at nadir corresponding to scene center
meta.data.clock_angle=str2double(fread(fid,8,'uint8=>char')); % Sensor clock angle as measured relative top sensor platform flight direction
meta.data.incidence=str2double(fread(fid,8,'uint8=>char')); % Incidence angle at scene center
fseek(fid,8,'cof'); % Blanks
meta.data.wavelength=str2double(fread(fid,16,'uint8=>char')); % Radar wavelength [m] = Nominal value
meta.data.mocomp=fread(fid,2,'uint8=>char').'; % Motion compensation indicator
meta.data.range_pulse_code=fread(fid,16,'uint8=>char').'; % Range pulse code
meta.data.range_pulse_amp_coef1=fread(fid,16,'uint8=>char').'; % Range pulse almplitude coefficient #1
meta.data.range_pulse_amp_coef2=fread(fid,16,'uint8=>char').'; % Range pulse almplitude coefficient #2
meta.data.range_pulse_amp_coef3=fread(fid,16,'uint8=>char').'; % Range pulse almplitude coefficient #3
meta.data.range_pulse_amp_coef4=fread(fid,16,'uint8=>char').'; % Range pulse almplitude coefficient #4
meta.data.range_pulse_amp_coef5=fread(fid,16,'uint8=>char').'; % Range pulse almplitude coefficient #5
meta.data.range_pulse_phs_coef1=fread(fid,16,'uint8=>char').'; % Range pulse phase coefficient #1
meta.data.range_pulse_phs_coef2=fread(fid,16,'uint8=>char').'; % Range pulse phase coefficient #2
meta.data.range_pulse_phs_coef3=fread(fid,16,'uint8=>char').'; % Range pulse phase coefficient #3
meta.data.range_pulse_phs_coef4=fread(fid,16,'uint8=>char').'; % Range pulse phase coefficient #4
meta.data.range_pulse_phs_coef5=fread(fid,16,'uint8=>char').'; % Range pulse phase coefficient #5
meta.data.chirp_index=fread(fid,8,'uint8=>char').'; % Downlinked data chirp extraction index
fseek(fid,8,'cof'); % Blanks
meta.data.sampling_rate=str2double(fread(fid,16,'uint8=>char')); % Sampling rate [MHz]
meta.data.range_gate=str2double(fread(fid,16,'uint8=>char')); % Range gate (early edge (in time) at the start of the image) [usec]
meta.data.pulse_width=str2double(fread(fid,16,'uint8=>char')); % Range pulse width [usec]
meta.data.baseband_flg=fread(fid,4,'uint8=>char').'; % Baseband conversion flag
meta.data.range_compressed_flg=fread(fid,4,'uint8=>char').'; % Range compressed flag
meta.data.rec_gain_like_pol=str2double(fread(fid,16,'uint8=>char')); % Receiver gain for like polarized at early edge at the start of the image
meta.data.rec_gain_cross_pol=str2double(fread(fid,16,'uint8=>char')); % Receiver gain for cross polarized at early edge at the start of the image
meta.data.quant_bit=fread(fid,8,'uint8=>char').'; % Quantization in bits per channel
meta.data.quant_desc=fread(fid,12,'uint8=>char').'; % Quantized descriptor
meta.data.dc_bias_i=str2double(fread(fid,16,'uint8=>char')); % DC Bias for I-component
meta.data.dc_bias_q=str2double(fread(fid,16,'uint8=>char')); % DC Bias for Q-component
meta.data.gain_imbalance=str2double(fread(fid,16,'uint8=>char')); % Gain imbalance for I And Q
fseek(fid,32,'cof'); % Blanks
meta.data.elec_bores=str2double(fread(fid,16,'uint8=>char')); % Electronic boresight
meta.data.mech_bores=str2double(fread(fid,16,'uint8=>char')); % Mechanical boresight
meta.data.echo_tracker=fread(fid,4,'uint8=>char').'; % Echo tracker on/off
meta.data.prf=str2double(fread(fid,16,'uint8=>char')); % PRF [mHz]
meta.data.ant_beam_2way_el=str2double(fread(fid,16,'uint8=>char')); % Two-way antenna beam width [deg] (Elevation)
meta.data.ant_beam_2way_az=str2double(fread(fid,16,'uint8=>char')); % Two-way antenna beam width [deg] (Azimuth)
meta.data.sat_time=fread(fid,16,'uint8=>char').'; % Satellite encoded binary time code (Tref)
meta.data.sat_clock=fread(fid,32,'uint8=>char').'; % Satellite clock time: Standard ground time of error time information (Tgref)
meta.data.sat_clock_inc=fread(fid,16,'uint8=>char').'; % Satellite clock increment [nsec]: Error time information of calculation satellite counter cycle
meta.data.proc_fac=fread(fid,16,'uint8=>char').'; % Processing facility ID
meta.data.proc_sys=fread(fid,8,'uint8=>char').'; % Processing system ID
meta.data.proc_ver=fread(fid,8,'uint8=>char').'; % Processing version ID
meta.data.proc_fac_code=fread(fid,16,'uint8=>char').'; % Processing code of processing facility
meta.data.proc_lvl_code=fread(fid,16,'uint8=>char').'; % Processing level code
meta.data.prod_type=fread(fid,32,'uint8=>char').'; % Product type specifier
meta.data.proc_alg=fread(fid,32,'uint8=>char').'; % Processing algorithm ID
meta.data.num_looks_az=str2double(fread(fid,16,'uint8=>char')); % Number of looks in azimuth
meta.data.num_looks_rng=str2double(fread(fid,16,'uint8=>char')); % Number of looks in range
meta.data.bw_per_look_az=str2double(fread(fid,16,'uint8=>char')); % Bandwidth per look in azimuth [Hz]
meta.data.bw_per_look_rng=str2double(fread(fid,16,'uint8=>char')); % Bandwidth per look in range [Hz]
meta.data.bw_az=str2double(fread(fid,16,'uint8=>char')); % Bandwidth in azimuth [Hz]
meta.data.bw_rng=str2double(fread(fid,16,'uint8=>char')); % Bandwidth in range [kHz]
meta.data.wgt_az=fread(fid,32,'uint8=>char').'; % Weighting function in azimuth
meta.data.wgt_rng=fread(fid,32,'uint8=>char').'; % Weighting function in range
meta.data.data_input_src=fread(fid,16,'uint8=>char').'; % Data input source
meta.data.res_grnd_rng=fread(fid,16,'uint8=>char'); % Resolution in ground range [m]
meta.data.res_az=fread(fid,16,'uint8=>char'); % Resolution in azimuth [m]
meta.data.rad_bias=fread(fid,16,'uint8=>char'); % Radiometric parameter (bias)
meta.data.rad_gain=fread(fid,16,'uint8=>char'); % Radiometric parameter (gain)
meta.data.at_dop_const=str2double(fread(fid,16,'uint8=>char')); % Along track Doppler frequency (center) constant term at early edge of image [Hz]
meta.data.at_dop_lin=str2double(fread(fid,16,'uint8=>char')); % Along track Doppler frequency (center) linear term at early edge of image [Hz]
meta.data.at_dop_quad=str2double(fread(fid,16,'uint8=>char')); % Along track Doppler frequency (center) quadratic term at early edge of image [Hz]
fseek(fid,16,'cof'); % Blanks
meta.data.xt_dop_const=str2double(fread(fid,16,'uint8=>char')); % Cross track Doppler frequency (center) constant term at early edge of image [Hz]
meta.data.xt_dop_lin=str2double(fread(fid,16,'uint8=>char')); % Cross track Doppler frequency (center) linear term at early edge of image [Hz]
meta.data.xt_dop_quad=str2double(fread(fid,16,'uint8=>char')); % Cross track Doppler frequency (center) quadratic term at early edge of image [Hz]
meta.data.time_dir_pixel=fread(fid,8,'uint8=>char').'; % Time direction indicator along pixel direction
meta.data.time_dir_line=fread(fid,8,'uint8=>char').'; % Time direction indicator along line direction
meta.data.at_dop_rate_const=str2double(fread(fid,16,'uint8=>char')); % Along track Doppler frequency rate constant term at early edge of image [Hz]
meta.data.at_dop_rate_lin=str2double(fread(fid,16,'uint8=>char')); % Along track Doppler frequency rate linear term at early edge of image [Hz]
meta.data.at_dop_rate_quad=str2double(fread(fid,16,'uint8=>char')); % Along track Doppler frequency rate quadratic term at early edge of image [Hz]
fseek(fid,16,'cof'); % Blanks
meta.data.xt_dop_rate_const=str2double(fread(fid,16,'uint8=>char')); % Cross track Doppler frequency rate constant term at early edge of image [Hz]
meta.data.xt_dop_rate_lin=str2double(fread(fid,16,'uint8=>char')); % Cross track Doppler frequency rate linear term at early edge of image [Hz]
meta.data.xt_dop_rate_quad=str2double(fread(fid,16,'uint8=>char')); % Cross track Doppler frequency rate quadratic term at early edge of image [Hz]
fseek(fid,16,'cof'); % Blanks
meta.data.line_constant=fread(fid,8,'uint8=>char').'; % Line constant indicator
meta.data.clutter_lock_flg=fread(fid,4,'uint8=>char').'; % Clutter lock applied flag
meta.data.autofocus_flg=fread(fid,4,'uint8=>char').'; % Auto-focusing applied flag
meta.data.line_spacing=str2double(fread(fid,16,'uint8=>char')); % Line spacing [m]
meta.data.pixel_spacing=str2double(fread(fid,16,'uint8=>char')); % Pixel spacing [m]
meta.data.rng_comp_des=fread(fid,16,'uint8=>char').'; % Processor range compression designator
meta.data.dop_freq_const=str2double(fread(fid,16,'uint8=>char')); % Doppler center frequency approximately coefficient constant term
meta.data.dop_freq_lin=str2double(fread(fid,16,'uint8=>char')); % Doppler center frequency approximately linear coefficient term
meta.data.cal_mode_loc_flag=fread(fid,4,'uint8=>char').'; % Calibration mode data location flag
meta.data.start_line_cal_start=fread(fid,8,'uint8=>char').'; % Start line number of calibration at the side of start
meta.data.end_line_cal_start=fread(fid,8,'uint8=>char').'; % End line number of calibration at the side of start
meta.data.start_line_cal_end=fread(fid,8,'uint8=>char').'; % Start line number of calibration at the side of end
meta.data.end_line_cal_end=fread(fid,8,'uint8=>char').'; % End line number of calibration at the side of end
meta.data.prf_switch=fread(fid,4,'uint8=>char').'; % PRF switching indicator
meta.data.prf_switch_line=fread(fid,8,'uint8=>char').'; % Line number of PRF switching
meta.data.beam_ctr_dir=str2double(fread(fid,16,'uint8=>char')); % The direction of a beam center in a scene center
meta.data.yaw_steer_flg=fread(fid,4,'uint8=>char').'; % Yaw steering mode flag
meta.data.param_table_num=fread(fid,4,'uint8=>char').'; % Parameter table number of automatically setting
meta.data.off_nadir=str2double(fread(fid,16,'uint8=>char')); % Nominal off nadir angle
meta.data.ant_beam_num=fread(fid,4,'uint8=>char').'; % Antenna beam number
fseek(fid,28,'cof'); % Blanks
meta.data.a0=str2double(fread(fid,20,'uint8=>char')); % Incidence angle constant term (a0)
meta.data.a1=str2double(fread(fid,20,'uint8=>char')); % Incidence angle linear coefficient term (a1)
meta.data.a2=str2double(fread(fid,20,'uint8=>char')); % Incidence angle quadratic coefficient term (a2)
meta.data.a3=str2double(fread(fid,20,'uint8=>char')); % Incidence angle cubic coefficient term (a3)
meta.data.a4=str2double(fread(fid,20,'uint8=>char')); % Incidence angle fourth coefficient term (a4)
meta.data.a5=str2double(fread(fid,20,'uint8=>char')); % Incidence angle fifth coefficient term (a5)
meta.data.num_annot=str2double(fread(fid,8,'uint8=>char')); % Number of annotation points (should be zero)
fseek(fid,8+64*32+26,'cof'); % Blanks

% Read map projection
if meta.num_map_rec % Not use for level 1.1 data
    fseek(fid,meta.map_len,'cof');
    % meta.map.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
    % meta.map.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
    % meta.map.rec_type=fread(fid,1,'uint8'); % Record type
    % meta.map.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
    % meta.map.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
    % meta.map.rec_length=fread(fid,1,'uint32'); % Record length
    % fseek(fid,16,'cof'); % Blanks
    % meta.map.proj=fread(fid,32,'uint8=>char').'; % Map projection
    % ...
end

% Read platform position data
pos_pos = ftell(fid);
meta.pos.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
meta.pos.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
meta.pos.rec_type=fread(fid,1,'uint8'); % Record type
meta.pos.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
meta.pos.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
meta.pos.rec_length=fread(fid,1,'uint32'); % Record length
meta.pos.orb_elem=fread(fid,32,'uint8=>char').'; % Orbital elements designator
meta.pos.pos_x=str2double(fread(fid,16,'uint8=>char')); % Position (x) [m]
meta.pos.pos_y=str2double(fread(fid,16,'uint8=>char')); % Position (y) [m]
meta.pos.pos_z=str2double(fread(fid,16,'uint8=>char')); % Position (z) [m]
meta.pos.vel_x=str2double(fread(fid,16,'uint8=>char')); % Velocity (x) [m]
meta.pos.vel_y=str2double(fread(fid,16,'uint8=>char')); % Velocity (y) [m]
meta.pos.vel_z=str2double(fread(fid,16,'uint8=>char')); % Velocity (z) [m]
meta.pos.num_pts=str2double(fread(fid,4,'uint8=>char')); % Number of data points (28)
meta.pos.year=str2double(fread(fid,4,'uint8=>char')).'; % Year of 1st point
meta.pos.month=str2double(fread(fid,4,'uint8=>char')).'; % Month of 1st point
meta.pos.day=str2double(fread(fid,4,'uint8=>char')).'; % Day of 1st point
meta.pos.day_in_year=str2double(fread(fid,4,'uint8=>char')).'; % Day in year of 1st point
meta.pos.sec=str2double(fread(fid,22,'uint8=>char')); % Seconds of day of 1st point
meta.pos.int=str2double(fread(fid,22,'uint8=>char')); % Time interval between data points
meta.pos.ref_coord_sys=fread(fid,64,'uint8=>char').'; % Reference coordinate system
meta.pos.greenwich_mean_hr_ang=fread(fid,22,'uint8=>char').'; % Greenwich mean hour angle
meta.pos.at_pos_err=str2double(fread(fid,16,'uint8=>char')); % Along track position error
meta.pos.ct_pos_err=str2double(fread(fid,16,'uint8=>char')); % Across track position error
meta.pos.rad_pos_err=str2double(fread(fid,16,'uint8=>char')); % Radial position error
meta.pos.at_vel_err=str2double(fread(fid,16,'uint8=>char')); % Along track velocity error
meta.pos.ct_vel_err=str2double(fread(fid,16,'uint8=>char')); % Across track velocity error
meta.pos.rad_vel_err=str2double(fread(fid,16,'uint8=>char')); % Radial velocity error
for i = 1:meta.pos.num_pts
    meta.pos.pts(i).pos_x=str2double(fread(fid,22,'uint8=>char')); % Position (x) [m]
    meta.pos.pts(i).pos_y=str2double(fread(fid,22,'uint8=>char')); % Position (y) [m]
    meta.pos.pts(i).pos_z=str2double(fread(fid,22,'uint8=>char')); % Position (z) [m]
    meta.pos.pts(i).vel_x=str2double(fread(fid,22,'uint8=>char')); % Velocity (x) [m]
    meta.pos.pts(i).vel_y=str2double(fread(fid,22,'uint8=>char')); % Velocity (y) [m]
    meta.pos.pts(i).vel_z=str2double(fread(fid,22,'uint8=>char')); % Velocity (z) [m]
end
fseek(fid,18,'cof'); % Blanks
meta.pos.leap_sec=fread(fid,1,'uint8=>char').'; % Occurrence flag of a leap second
fseek(fid,pos_pos+meta.pos.rec_length,'bof'); % Blanks

% Read attitude data
att_pos = ftell(fid);
meta.att.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
meta.att.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
meta.att.rec_type=fread(fid,1,'uint8'); % Record type
meta.att.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
meta.att.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
meta.att.rec_length=fread(fid,1,'uint32'); % Record length
meta.att.num_pts=str2double(fread(fid,4,'uint8=>char')); % Number of points
for i = 1:meta.att.num_pts
    meta.att.pts(i).day_year=str2double(fread(fid,4,'uint8=>char')); % Day of the year
    meta.att.pts(i).msec_day=str2double(fread(fid,8,'uint8=>char')); % Millisecond of the day
    meta.att.pts(i).pitch_flag=fread(fid,4,'uint8=>char').'; % Pitch data quality flag
    meta.att.pts(i).roll_flag=fread(fid,4,'uint8=>char').'; % Roll data quality flag
    meta.att.pts(i).yaw_flag=fread(fid,4,'uint8=>char').'; % Yaw data quality flag
    meta.att.pts(i).pitch=str2double(fread(fid,14,'uint8=>char')); % Pitch [deg]
    meta.att.pts(i).roll=str2double(fread(fid,14,'uint8=>char')); % Roll [deg]
    meta.att.pts(i).yaw=str2double(fread(fid,14,'uint8=>char')); % Yaw [deg]
    meta.att.pts(i).pitch_rate_flag=fread(fid,4,'uint8=>char').'; % Pitch rate data quality flag
    meta.att.pts(i).roll_rate_flag=fread(fid,4,'uint8=>char').'; % Roll rate data quality flag
    meta.att.pts(i).yaw_rate_flag=fread(fid,4,'uint8=>char').'; % Yaw rate data quality flag
    meta.att.pts(i).pitch_rate=str2double(fread(fid,14,'uint8=>char')); % Pitch rate
    meta.att.pts(i).roll_rate=str2double(fread(fid,14,'uint8=>char')); % Roll rate
    meta.att.pts(i).yaw_rate=str2double(fread(fid,14,'uint8=>char')); % Yaw rate
end
fseek(fid,att_pos+meta.att.rec_length,'bof'); % Blanks
%fseek(fid,meta.att.rec_length-(16+120*meta.att.num_pts),'cof'); % Blanks

% Read radiometric data
meta.rad.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
meta.rad.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
meta.rad.rec_type=fread(fid,1,'uint8'); % Record type
meta.rad.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
meta.rad.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
meta.rad.rec_length=fread(fid,1,'uint32'); % Record length
meta.rad.seq_num=str2double(fread(fid,4,'uint8=>char')); % Radiometric data records sequence number
meta.rad.num_pts=str2double(fread(fid,4,'uint8=>char')); % Number of radiometric fields
meta.rad.cal_factor=str2double(fread(fid,16,'uint8=>char')); % Calibration factor
meta.rad.dt = zeros(8,1);
for i = 1:8
    meta.rad.dt(i)=str2double(fread(fid,16,'uint8=>char')); % Transmission distortion matrix for high-sensitive/fine modes (full (quad) pol)
end
meta.rad.dt = reshape(meta.rad.dt(1:2:end) + 1j*meta.rad.dt(2:2:end),[2 2]).';
meta.rad.dr = zeros(8,1);
for i = 1:8
    meta.rad.dr(i)=str2double(fread(fid,16,'uint8=>char')); % Receive distortion matrix for high-sensitive/fine modes (full (quad) pol)
end
meta.rad.dr = reshape(meta.rad.dr(1:2:end) + 1j*meta.rad.dr(2:2:end),[2 2]).';
fseek(fid,9568,'cof'); % Blanks

% Read data quality records
data_qual_pos = ftell(fid);
meta.data_qual.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
meta.data_qual.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
meta.data_qual.rec_type=fread(fid,1,'uint8'); % Record type
meta.data_qual.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
meta.data_qual.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
meta.data_qual.rec_length=fread(fid,1,'uint32'); % Record length
meta.data_qual.dq_rec_num=fread(fid,4,'uint8=>char').'; % Radiometric data records sequence number
meta.data_qual.chan_id=fread(fid,4,'uint8=>char').'; % SAR Channel ID
meta.data_qual.date=fread(fid,6,'uint8=>char').'; % Date of last calibration update
meta.data_qual.num_chans=str2double(fread(fid,4,'uint8=>char')); % Number of channels
meta.data_qual.islr=str2double(fread(fid,16,'uint8=>char')); % ISLR
meta.data_qual.pslr=str2double(fread(fid,16,'uint8=>char')); % PSLR
meta.data_qual.aar=str2double(fread(fid,16,'uint8=>char')); % Azimuth ambiguity rate (nominal value)
meta.data_qual.rar=str2double(fread(fid,16,'uint8=>char')); % Range ambiguity rate (nominal value)
meta.data_qual.snr=str2double(fread(fid,16,'uint8=>char')); % Estimate of SNR [dB]
meta.data_qual.ber=fread(fid,16,'uint8=>char').'; % BER (actual value)
meta.data_qual.sr_res=str2double(fread(fid,16,'uint8=>char')); % Slant range resolution (nominal value) [m]
meta.data_qual.az_res=str2double(fread(fid,16,'uint8=>char')); % Azimuth resolution (nominal value) [m]
meta.data_qual.rad_res=fread(fid,16,'uint8=>char').'; % Radiometric resolution [dB]
meta.data_qual.dyn_rng=fread(fid,16,'uint8=>char').'; % Instantateous dynamic range [dB]
meta.data_qual.abs_cal_mag=fread(fid,16,'uint8=>char').'; % Nominal absolute radiometric calibration magnitude uncertainty
meta.data_qual.abs_cal_phs=fread(fid,16,'uint8=>char').'; % Nominal absolute radiometric calibration phase uncertainty
for i = 1:meta.data_qual.num_chans
    meta.data_qual.rel_cal_mag(i)=str2double(fread(fid,16,'uint8=>char')); % Nominal relative radiometric calibration magnitude uncertainty
    meta.data_qual.rel_cal_phs(i)=str2double(fread(fid,16,'uint8=>char')); % Nominal relative radiometric calibration phase uncertainty
end
fseek(fid,512-(meta.data_qual.num_chans*32),'cof'); % Blanks
meta.data_qual.abs_err_at=fread(fid,16,'uint8=>char').'; % Absolute location error along track
meta.data_qual.abs_err_ct=fread(fid,16,'uint8=>char').'; % Absolute location error cross track
meta.data_qual.distort_line=fread(fid,16,'uint8=>char').'; % Geometric distortion scale in line direction
meta.data_qual.distort_pixel=fread(fid,16,'uint8=>char').'; % Geometric distortion scale in pixel direction
meta.data_qual.distort_skew=fread(fid,16,'uint8=>char').'; % Geometric distortion skew
meta.data_qual.orient_err=fread(fid,16,'uint8=>char').'; % Scene orientation error
for i = 1:meta.data_qual.num_chans
    meta.data_qual.at_misreg_err(i)=str2double(fread(fid,16,'uint8=>char')); % Along track relative misregistration error
    meta.data_qual.ct_misreg_err(i)=str2double(fread(fid,16,'uint8=>char')); % Cross track relative misregistration error
end
fseek(fid,data_qual_pos+meta.data_qual.rec_length,'bof'); % Blanks
%fseek(fid,534,'cof'); % Blanks

% Facility related records
for i=1:20
    fac_pos = ftell(fid);
    meta.fac{i}.rec_seq_num=fread(fid,1,'uint32'); % Record sequence number
    if feof(fid)
        meta.fac = meta.fac(1:(i-1));
        break;
    end
    meta.fac{i}.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
    meta.fac{i}.rec_type=fread(fid,1,'uint8'); % Record type
    meta.fac{i}.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
    meta.fac{i}.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
    meta.fac{i}.rec_length=fread(fid,1,'uint32'); % Record length
    meta.fac{i}.fac_seq_num=str2double(fread(fid,4,'uint8=>char')); % Facility related data record sequence number
    if meta.fac{i}.rec_seq_num == 11
        for j = 1:10
            meta.fac{i}.latlon2pix(j)=str2double(fread(fid,20,'uint8=>char')); % Coefficients to convert from map project (E, N) to pixel
        end
        for j = 1:10
            meta.fac{i}.latlon2line(j)=str2double(fread(fid,20,'uint8=>char')); % Coefficients to convert from map project (E, N) to line
        end
        meta.fac{i}.cal_mode_data_loc_flg=fread(fid,4,'uint8=>char').'; % Calibration mode data location flag
        meta.fac{i}.start_line_upper=fread(fid,8,'uint8=>char').'; % Start line number of calibration at upper image
        meta.fac{i}.end_line_upper=fread(fid,8,'uint8=>char').'; % Start line number of calibration at upper image
        meta.fac{i}.start_line_bottom=fread(fid,8,'uint8=>char').'; % Start line number of calibration at upper image
        meta.fac{i}.end_line_bottom=fread(fid,8,'uint8=>char').'; % Start line number of calibration at upper image
        meta.fac{i}.prf_switch_flag=fread(fid,4,'uint8=>char').'; % PRF switching flag
        meta.fac{i}.prf_switch_line=fread(fid,8,'uint8=>char').'; % Start line number of PRF switching
        fseek(fid,8,'cof'); % Blanks
        meta.fac{i}.num_loss_lines_10=fread(fid,8,'uint8=>char').'; % Number of loss lines (Level 1.0)
        meta.fac{i}.num_loss_lines_11=fread(fid,8,'uint8=>char').'; % Number of loss lines (Level 1.1, 1.5, 3.1)
        fseek(fid,312,'cof'); % Blanks
        fseek(fid,224,'cof'); % System reserve
        for j = 1:25
            meta.fac{i}.a(j)=str2double(fread(fid,20,'uint8=>char')); % Coefficients to convert from pixel and line to latitude
        end
        for j = 1:25
            meta.fac{i}.b(j)=str2double(fread(fid,20,'uint8=>char')); % Coefficients to convert from pixel and line to longitude
        end
        meta.fac{i}.origin_pixel=str2double(fread(fid,20,'uint8=>char')); % Origin pixel
        meta.fac{i}.origin_line=str2double(fread(fid,20,'uint8=>char')); % Origin line
        for j = 1:25
            meta.fac{i}.c(j)=str2double(fread(fid,20,'uint8=>char')); % Coefficients to convert from latitude and longitude to pixel
        end
        for j = 1:25
            meta.fac{i}.d(j)=str2double(fread(fid,20,'uint8=>char')); % Coefficients to convert from latitude and longitude to line
        end
        meta.fac{i}.origin_lat=str2double(fread(fid,20,'uint8=>char')); % Origin latitude
        meta.fac{i}.origin_lon=str2double(fread(fid,20,'uint8=>char')); % Origin longitude
    else
        % Unclear how to parse dummy data, determined ephemeris, time error
        % information, or coordinate conversion information that are stored in
        % the facility related records.
    end
    fseek(fid,fac_pos+meta.fac{i}.rec_length,'bof'); % Blanks
end

% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////