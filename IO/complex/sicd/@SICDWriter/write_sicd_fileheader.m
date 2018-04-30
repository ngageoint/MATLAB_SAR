function write_sicd_fileheader(obj)
%WRITE_SICD_FILEHEADER Writes the NITF file header
%
% It's not expected that the typical user would call this function, rather
% it is called as part of a bigger NITF/SICD writer.
%
% Inputs:
%       obj:        This is the 'hidden' object-specific handle.  It is
%                   similar to 'this' in Java.
%
% References:
%       NGA.STND.0024-2_1.0, Sensor Independent Complex Data (SICD), Volume
%          2, File Format Description Document, 2011-08-01
%       MIL-STD-2500C, Department of Defense Interface Standard
%                      National Imagery Transmission Format version 2.1
%                      01 May 2006
%
% Written by: Tom Krauss and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Note, the hard-coded numbers in the writeN calls are NITF 2.1 field
% widths (in characters).

% Always writing to NITF version 2.1
fwriten(obj.FID, 'NITF02.10', 9);  % FHDR and FVER
% Complexity level depends on the "size" of the image in bytes.
image_data_size = sum(obj.NumRowsIS) * obj.BytesPerRow; 
if image_data_size < (50*(1024*1024)) % Less than 50MB
    complexity = 3;
elseif image_data_size < (1*(1024*1024*1024)) % Less than 1GB
    complexity = 5;
elseif image_data_size < (2*(1024*1024*1024)) % Less than 2GB
    complexity = 6;
else
    complexity = 7;
end
fwriten(obj.FID, uint8(complexity),  2);  % CLEVEL
fwriten(obj.FID, 'BF01', 4);  % STYPE
fwriten(obj.FID, 'Unknown', 10); % OSTAID (not supposed to be blank)

if isfield(obj.sicdmeta,'ImageCreation')&&isfield(obj.sicdmeta.ImageCreation,'DateTime')
    fdt = datestr(obj.sicdmeta.ImageCreation.DateTime,'yyyymmddHHMMSS'); % Creation time of original image
else
    fdt = datestr(now,'yyyymmddHHMMSS'); % Creation time of this NITF
end
fwriten(obj.FID, fdt, 14); % FDT

if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'CoreName')
    FTITLE = obj.sicdmeta.CollectionInfo.CoreName;
else
    FTITLE = 'Unknown';
end
fwriten(obj.FID, ['SICD: ' FTITLE], 80); % FTITLE

obj.write_sicd_security_tags(); % Security tags go here
fwriten(obj.FID, '00000', 5);  % FSCOP - File Copy Number
fwriten(obj.FID, '00000', 5);  % FSCPYS - File Number of Copies
fwriten(obj.FID, '0', 1); % ENCRYP - Encryption. (0=none)
fwrite(obj.FID, 0, 'uint8'); % FBKGC 3-byte binary write for
fwrite(obj.FID, 0, 'uint8'); % background color
fwrite(obj.FID, 0, 'uint8'); % (red, green, blue)
fwriten(obj.FID, '', 24); % ONAME
fwriten(obj.FID, '', 18); % OPHONE

% File and header lengths
fwriten(obj.FID, uint64(obj.NITF_header_length + ...
    (obj.IS_SUBHEADER_LENGTH * obj.NumIS) + ...
    image_data_size + ...
    obj.DES_HEADER_LENGTH) + ...
    length(obj.DES_data), 12);
fwriten(obj.FID, uint64(obj.NITF_header_length), 6);

% Image Segment description
fwriten(obj.FID, uint16(obj.NumIS), 3);
for i = 1:obj.NumIS
    fwriten(obj.FID, uint64(obj.IS_SUBHEADER_LENGTH), 6);
    fwriten(obj.FID, uint64(obj.NumRowsIS(i) * obj.BytesPerRow), 10);
end

% Graphic segments (not allowed in SICD)
fwriten(obj.FID, uint64(0), 3);

% Reserved extension segments (not allowed in SICD)
fwriten(obj.FID, uint64(0), 3);

% Text Segments - (Not generally used in SICD)
fwriten(obj.FID, uint64(0), 3);

% Data Extension segments.  This is where the SICD XML goes.  We'll make
% the XML string here so we can get its length even though we won't write
% the string until later.
fwriten(obj.FID, uint64(1), 3);
fwriten(obj.FID, uint64(obj.DES_HEADER_LENGTH), 4);
fwriten(obj.FID, uint64(length(obj.DES_data)), 9);   % Length of DES data

% Reserved Extension segments (not allowed in SICD)
fwriten(obj.FID, uint64(0), 3);

% User defined headers (Not generally used in SICD)
fwriten(obj.FID, uint64(0), 5);

% Extended headers (Not generally used in SICD)
fwriten(obj.FID, uint64(0), 5);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////