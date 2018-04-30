function [ thumb_out ] = read_sicd01_thumb( filename, meta_in )
%READ_SICD01_THUMB Read thumbnail image from Sensor Independent Complex
%Data (SICD) file, version 0.1
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Use metadata is already read in, use current copy.  Otherwise, read it in again.
if(nargin<2)
    meta_in=read_sicd01_meta( filename );
end

% Skip to appropriate place in file to read JPEG thumbnail
% non-xml header + xml + data + 6 bytes (3 blank lines to seperate header,xml,data,jpeg)
fid=fopen(filename,'r','b');
fseek(fid,meta_in.header_length_in_bytes+meta_in.header.xml_content_length+meta_in.header.data_content_length+6,'bof');
jpegdata=fread(fid,meta_in.header.additional_data_content_length,'uchar',0)';
fclose(fid);
keepname=[tempname '.jpg']; % Store data in temporary file and use Matlab's imread to decode JPEG
fid2=fopen(keepname,'wb');
if(fwrite(fid2,jpegdata,'uchar')==length(jpegdata))
    fclose(fid2);
    try
        thumb_out=imread(keepname);
    catch
        warning('READ_SICD01_THUMB:InvalidJPEG','Invalid JPEG thumbnail data contained in file.');
        thumb_out=[];
    end
else
    fclose(fid2);
    warning('READ_SICD01_THUMB:TempFileError','Unable to write temporary file necessary to decode JPEG thumbnail.');
    thumb_out=[];
end
delete(keepname);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////