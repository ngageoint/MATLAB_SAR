function [ boolout ] = ispalsar2( filename )
%ISPALSAR2 ALOS PALSAR 2 file format
%
% Checks to see if file is any of the files in the volume, leader, image,
% trailer package.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

try
    % This will parse the first few bytes for the header of any of the file
    % types in the package.
    fid=fopen(filename,'r','b');
    meta.rec_num=fread(fid,1,'uint32'); % Record sequence number
    meta.rec_subtype1=fread(fid,1,'uint8'); % Record subtype code
    meta.rec_type=fread(fid,1,'uint8'); % Record type
    meta.rec_subtype2=fread(fid,1,'uint8'); % Record subtype code
    meta.rec_subtype3=fread(fid,1,'uint8'); % Record subtype code
    meta.rec_length=fread(fid,1,'uint32'); % Record length
    meta.ascii_ebcdic=fread(fid,2,'uint8=>char').'; % Record length
    fseek(fid,2,'cof'); % Blanks
    meta.doc_id=fread(fid,12,'uint8=>char').'; % Superstructure format control document ID
    fclose(fid);

    % This only strictly and explicitly checks for CEOS format, not
    % necessarily for ALOS-2.  For instance, ALOS-1 and RADARSAT-1 are also
    % a SARs that deliver their data in CEOS format (but these both use a
    % different doc_id string, so hopefully this is OK).
    boolout = meta.rec_num==1 && ...
        any(meta.rec_subtype1==[192, 11, 50, 63]) && ...
        meta.rec_type==192 && ...
        meta.rec_subtype2==18 && ...
        meta.rec_subtype3==18 && ...
        any(meta.rec_length==[360, 720]) && ...
        strncmp(meta.doc_id, 'CEOS-SAR', 8);
catch
    boolout = false;
end
% If a non-image file (VOL/LED/TRL) file is passed, make sure an image file
% also exist.
if boolout && (meta.rec_subtype1 ~= 50)
    filenames = palsar2_files(filename, false);
    boolout = isfield(filenames, 'IMG');
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////