function [ meta ] = read_ceos_trl_meta( filename )
%READ_CEOS_TRL_META Read CEOS SAR trailer file
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open file
fid=fopen(filename,'r','b');

% Read volume descriptor records
meta.record_sequency_number=fread(fid,1,'uint32'); % Record sequence number
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
% The following should all be zero
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
% Low resolution records
meta.num_low_res_rec=str2double(fread(fid,6,'uint8=>char')); % Number of low resolution records
for i = 1:meta.num_low_res_rec
    meta.low_res(i).len=str2double(fread(fid,8,'uint8=>char')); % Low resolution image data record length
    meta.low_res(i).pixels=str2double(fread(fid,6,'uint8=>char')); % Number of pixels of low resolution image data
    meta.low_res(i).lines=str2double(fread(fid,6,'uint8=>char')); % Number of lines of low resolution image data
    meta.low_res(i).bytes=str2double(fread(fid,6,'uint8=>char')); % Number of bytes per one sample of low resolution image data
end
fseek(fid,720,'bof'); % Blanks
% Seems to be an array the size of the low resolution image on the end of
% this file, but it doesn't seem to contain any data.
% meta.low_res(i).data = fread(fid,[meta.low_res(i).pixels, meta.low_res(i).lines],'uint16');

% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////