function [ meta ] = read_ceos_vol_meta( filename )
%READ_CEOS_VOL_META Read CEOS volume directory file
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
meta.phys_vol_id=fread(fid,16,'uint8=>char').'; % Physical volume ID
meta.log_vol_id=fread(fid,16,'uint8=>char').'; % Logical volume ID
meta.vol_set_id=fread(fid,16,'uint8=>char').'; % Volume set ID
meta.num_phys_vol=fread(fid,2,'uint8=>char').'; % Total number of physical volumes in the logical volume
meta.phys_seq_first=fread(fid,2,'uint8=>char').'; % Physical volume sequence number of the first tape
meta.phys_seq_last=fread(fid,2,'uint8=>char').'; % Physical volume sequence number of the last tape
meta.phys_seq_cur=fread(fid,2,'uint8=>char').'; % Physical volume sequence number of the current tape
meta.filen_num=fread(fid,4,'uint8=>char').'; % File number in the logical volume follows volume directory file
meta.log_vol=fread(fid,4,'uint8=>char').'; % Logical volume within a volume set
meta.log_vol_num=fread(fid,4,'uint8=>char').'; % Logical volume number with physical volume
meta.log_vol_create_date=fread(fid,8,'uint8=>char').'; % Logical volume creation date
meta.log_vol_create_time=fread(fid,8,'uint8=>char').'; % Logical volume creation time
meta.log_vol_co=fread(fid,12,'uint8=>char').'; % Logical volume generation country
meta.log_vol_agency=fread(fid,8,'uint8=>char').'; % Logical volume generation agency
meta.log_vol_facility=fread(fid,12,'uint8=>char').'; % Logical volume generation facility
meta.num_file_ptr=str2double(fread(fid,4,'uint8=>char')); % Number of file pointer records in volume directory
meta.num_text_rec=str2double(fread(fid,4,'uint8=>char')); % Number of text records in volume directory
fseek(fid,192,'cof'); % Blanks

% Read file pointer records
for i = 1:meta.num_file_ptr
    meta.file(i).rec_num = fread(fid,1,'uint32'); % Record number
    meta.file(i).rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
    meta.file(i).rec_type=fread(fid,1,'uint8'); % Record type
    meta.file(i).rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
    meta.file(i).rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
    meta.file(i).rec_length=fread(fid,1,'uint32'); % Record length
    meta.file(i).ascii_ebcdic=fread(fid,2,'uint8=>char').'; % Record length
    fseek(fid,2,'cof'); % Blanks
    meta.file(i).num=fread(fid,4,'uint8=>char').'; % Referenced file number
    meta.file(i).name=fread(fid,16,'uint8=>char').'; % Referenced file name ID
    meta.file(i).class=fread(fid,28,'uint8=>char').'; % Referenced file class
    meta.file(i).class_code=fread(fid,4,'uint8=>char').'; % Referenced file class code
    meta.file(i).type=fread(fid,28,'uint8=>char').'; % Referenced file type
    meta.file(i).type_code=fread(fid,4,'uint8=>char').'; % Referenced file type code
    meta.file(i).num_recs=str2double(fread(fid,8,'uint8=>char')); % Number of records in referenced file
    meta.file(i).len_first_rec=str2double(fread(fid,8,'uint8=>char')); % Length of first record in referenced file
    meta.file(i).max_rec_len=str2double(fread(fid,8,'uint8=>char')); % Maximum record length in referenced file
    meta.file(i).rec_len_type=fread(fid,12,'uint8=>char').'; % Refrenced file record length type
    meta.file(i).rec_len_type_code=fread(fid,4,'uint8=>char').'; % Refrenced file record length type code
    meta.file(i).phys_vol_first=str2double(fread(fid,2,'uint8=>char')); % The number of the physical volue set containing the first record of the file
    meta.file(i).phys_vol_last=str2double(fread(fid,2,'uint8=>char')); % The number of the physical volue set containing the last record of the file
    meta.file(i).rec_num_first=str2double(fread(fid,8,'uint8=>char')); % Record number of the first record appearing on this physical volume
    meta.file(i).rec_num_last=str2double(fread(fid,8,'uint8=>char')); % Record number of the last record appearing on this physical volume
    fseek(fid,200,'cof'); % Blanks
end

% Read text records
for i = 1:meta.num_text_rec
    meta.text(i).rec_num = fread(fid,1,'uint32'); % Record number
    meta.text(i).rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
    meta.text(i).rec_type=fread(fid,1,'uint8'); % Record type
    meta.text(i).rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
    meta.text(i).rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
    meta.text(i).rec_length=fread(fid,1,'uint32'); % Record length
    meta.text(i).ascii_ebcdic=fread(fid,2,'uint8=>char').'; % Record length
    fseek(fid,2,'cof'); % Blanks
    meta.text(i).prod_id=fread(fid,40,'uint8=>char').'; % Product ID
    meta.text(i).location=fread(fid,60,'uint8=>char').'; % Location and date/time of product creation
    meta.text(i).phys_id=fread(fid,40,'uint8=>char').'; % Physical tape ID
    meta.text(i).scene_id=fread(fid,40,'uint8=>char').'; % Scene ID
    meta.text(i).scene_loc_id=fread(fid,40,'uint8=>char').'; % Scene location ID
    fseek(fid,124,'cof');
end

% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////