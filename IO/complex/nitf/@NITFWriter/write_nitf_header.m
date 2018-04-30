function write_nitf_header(obj)
%WRITE_NITF_HEADER Writes the NITF file header
%
% It's not expected that the typical user would call this function, rather
% it is called as part of a bigger NITF/SICD writer.
%
% Inputs:
%       obj:        This is the 'hidden' object-specific handle.  It is
%                   similar to 'this' in Java.
%
% References:
%       MIL-STD-2500C, Department of Defense Interface Standard
%                      National Imagery Transmission Format version 2.1
%                      01 May 2006
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%
% Note, the hard-coded numbers in the writeN calls are NITF field widths
% (in characters).  If the spec changes or we choose to use something other
% than NITF 2.1 these values must change.
%

% Always writing to NITF version 2.1
fwriten(obj.FID, 'NITF02.10', 9);  % FHDR and FVER
% Complexity level depends on the "size" of the image.
% Size of image data.  Note, we're converting the integer row/col counts to
% double because the multiplication may overflow normal Matlab integer
% arithmetic.  We also need to get the pixel size which is
if obj.image_data_size < (50*(1024*1024)) % Less than 50MB
    complexity = 3;
elseif obj.image_data_size < (1*(1024*1024*1024)) % Less than 1GB
    complexity = 5;
elseif obj.image_data_size < (2*(1024*1024*1024)) % Less than 2GB
    complexity = 6;
else
    complexity = 7;
end
fwriten(obj.FID, uint8(complexity),  2);  % CLEVEL
fwriten(obj.FID, 'BF01',      4);  % STYPE
fwriten(obj.FID, 'Unknown',   10); % OSTAID (not supposed to be blank)

if isfield(obj.sicdmeta,'ImageCreation')&&isfield(obj.sicdmeta.ImageCreation,'DateTime')
    fdt = datestr(obj.sicdmeta.ImageCreation.DateTime,'yyyymmddHHMMSS'); % Creation time of original image
else
    fdt = datestr(now,'yyyymmddHHMMSS'); % Creation time of this NITF
end
fwriten(obj.FID, fdt,         14); % FDT

if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'CoreName')
    FTITLE = obj.sicdmeta.CollectionInfo.CoreName;
else
    FTITLE = 'Unknown';
end
if obj.is_complex
    FTITLE = ['SICD: ' FTITLE];
end
fwriten(obj.FID, FTITLE, 80); % FTITLE

obj.write_nitf_security_tags(); % Security tags go here
fwriten(obj.FID,   '00000', 5);  % FSCOP - File Copy Number
fwriten(obj.FID,   '00000', 5);  % FSCPYS - File Number of Copies
fwriten(obj.FID,   '0', 1); % ENCRYP - Encryption. (0=none)
fwrite(obj.FID, 0, 'uint8'); % FBKGC 3-byte binary write for
fwrite(obj.FID, 0, 'uint8'); % background color
fwrite(obj.FID, 0, 'uint8'); % (red, green, blue)
fwriten(obj.FID,   '', 24); % ONAME
fwriten(obj.FID,   '', 18); % OPHONE

% File and header lengths

fwriten(obj.FID, uint64(obj.total_file_length), 12);   % File length
fwriten(obj.FID, uint64(obj.NITF_header_length), 6);   % Header length - fixed size for NITF 2.1

% Image Segment description
fwriten(obj.FID, uint64(1), 3);    % Always going to write a single image segment.
fwriten(obj.FID, uint64(obj.IS_subheader_size), 6);  % Length of image subheader - fixed size (512) for NITF 2.1
fwriten(obj.FID, uint64(obj.image_data_size), 10);

% Graphic segments (none)
fwriten(obj.FID, uint64(0), 3);

% Reserved extension segments (none)
fwriten(obj.FID, uint64(0), 3);

% Text Segments (none)
fwriten(obj.FID, uint64(0), 3);

% Data Extension segments.  This is where the SICD XML goes.  We'll make
% the XML string here so we can get its length even though we won't write
% the string until later.
fwriten(obj.FID, uint64(1), 3);
fwriten(obj.FID, uint64(obj.DES_header_length), 4);
fwriten(obj.FID, uint64(length(obj.DES_data)), 9);   % Length of DES data

% Reserved Extension segments (none)
fwriten(obj.FID, uint64(0), 3);

% User defined headers (none)
fwriten(obj.FID, uint64(0), 5);

% Extended headers (none)
fwriten(obj.FID, uint64(0), 5);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////