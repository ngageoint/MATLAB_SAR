function lengthWritten = write_nitf_imsubhdr(obj)
%WRITE_NITF_IMSUBHDR Writes the NITF image segment sub-header
%
% It's not expected that the typical user would call this function, rather
% it is called as part of a bigger NITF/SICD writer.
%
% Inputs:
%       obj:       This is the 'hidden' object-specific handle.  It is
%                  similar to 'this' in Java.
%
% Outputs:
%       lengthWritten: The number of bytes written to disk
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


% Keep track of the starting position.  We'll use this to compute the
% length (in bytes) written.
fileStartPosition = ftell(obj.FID);

% Always writing version 2.1 NITF

fwriten(obj.FID, 'IM',      2);
if obj.is_complex
    IID1 = 'SICD000';
else
    IID1 = 'Unknown';
end
fwriten(obj.FID, IID1, 10); % IID1

if isfield(obj.sicdmeta,'ImageCreation')&&isfield(obj.sicdmeta.ImageCreation,'DateTime')
    fdt = datestr(obj.sicdmeta.ImageCreation.DateTime,'yyyymmddHHMMSS');
else
    fdt = '';
end
fwriten(obj.FID, fdt,       14); % IDATIM

fwriten(obj.FID, '',        17); % TGTID

if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'CoreName')
    IID2 = obj.sicdmeta.CollectionInfo.CoreName;
else
    IID2 = 'Unknown';
end
if obj.is_complex
    IID2 = ['SICD: ' IID2];
end
fwriten(obj.FID, IID2,        80); % IID2

obj.write_nitf_security_tags();  % Security tags go here
fwriten(obj.FID, '0',       1);  % ENCRYPT
if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'CollectorName')
    ISORCE = obj.sicdmeta.CollectionInfo.CollectorName;
else
    ISORCE = '';
end
fwriten(obj.FID, ISORCE,      42); % ISORCE
fwriten(obj.FID, uint64(obj.sicdmeta.ImageData.NumRows),  8); % NROWS
fwriten(obj.FID, uint64(obj.sicdmeta.ImageData.NumCols),  8); % NCOLS

if strncmp(obj.data_type,'int',3)
  fwriten(obj.FID, 'SI', 3);       % PVTYPE (here signed integer)
elseif strncmp(obj.data_type,'uint',4)
  fwriten(obj.FID, 'INT', 3);      % PVTYPE (here unsigned integer)
elseif strncmp(obj.data_type,'float',5)
  fwriten(obj.FID, 'R', 3);        % PVTYPE (here real)
end
if obj.is_complex
    fwriten(obj.FID, 'NODISPLY',8);  % IREP
elseif ~isempty(obj.LUT)
    fwriten(obj.FID, 'RGB/LUT', 8);  % IREP
else
    fwriten(obj.FID, 'MONO',    8);  % IREP
end
fwriten(obj.FID, 'SAR',     8);  % ICAT
element_size=obj.data_size*8; % Element size in bits
fwriten(obj.FID, uint64(element_size),2);  % ABPP
fwriten(obj.FID, 'R',       1);  % PJUST
fwriten(obj.FID, 'G',       1);  % ICORDS
if isfield(obj.sicdmeta,'GeoData')&&isfield(obj.sicdmeta.GeoData,'ImageCorners')
    geo_format = {'num_units',3,'include_symbols',false};
    geolo = [latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lat,'lat',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lon,'lon',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lat,'lat',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lon,'lon',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lat,'lat',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lon,'lon',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lat,'lat',geo_format{:}) ...
        latlonstr(obj.sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lon,'lon',geo_format{:})];
else
    geolo = '';
end
fwriten(obj.FID, geolo,     60); % IGEOLO

fwriten(obj.FID, uint64(0), 1);  % NICOM (# image comments)
fwriten(obj.FID, 'NC',      2);  % IC (No image compression)
fwriten(obj.FID, uint64(obj.is_complex+1), 1);  % NBANDS

if obj.is_complex
    % 'Band' 1: in phase
    fwriten(obj.FID, '',  2); % IREPBAND1
    fwriten(obj.FID, 'I', 6); % ISUBCAT1
    fwriten(obj.FID, 'N', 1); % IFC1
    fwriten(obj.FID, '',  3); % IMFLT1
    fwriten(obj.FID, '0', 1); % NLUTS1

    % 'Band' 2: quadrature
    fwriten(obj.FID, '',  2); % IREPBAND2
    fwriten(obj.FID, 'Q', 6); % ISUBCAT2
    fwriten(obj.FID, 'N', 1); % IFC2
    fwriten(obj.FID, '',  3); % IMFLT2
    fwriten(obj.FID, '0', 1); % NLUTS2
elseif ~isempty(obj.LUT)
    fwriten(obj.FID, 'LU',    2); % IREPBAND1
    fwriten(obj.FID, '',      6); % ISUBCAT1
    fwriten(obj.FID, 'N',     1); % IFC1
    fwriten(obj.FID, '',      3); % IMFLT1
    fwriten(obj.FID, '3',     1); % NLUTS1
    fwriten(obj.FID, '00256', 5); % NELUT1
    fwrite(obj.FID, round(255*obj.LUT), 'uint8'); % LUTDnm
else
    fwriten(obj.FID, 'M', 2);  % IREPBAND1 (here mono)
    fwriten(obj.FID, 'M', 6);  % ISUBCAT1 (here magnitude)
    fwriten(obj.FID, 'N', 1);  % IFC1
    fwriten(obj.FID, '',  3);  % IMFLT1
    fwriten(obj.FID, '0', 1);  % NLUTS1
end

fwriten(obj.FID, '0', 1); % ISYNC, always 0
fwriten(obj.FID, 'P', 1);  % IMODE (band interleaved by pixel)
fwriten(obj.FID, uint64(1), 4);  % NBPR - blocks per row
fwriten(obj.FID, uint64(1), 4);  % NBPC - blocks per column
if obj.sicdmeta.ImageData.NumCols>8192
    NPPBH = 0; % (zero means "use NCOLS")
else
    NPPBH = obj.sicdmeta.ImageData.NumCols;
end
fwriten(obj.FID, uint64(NPPBH), 4);  % NPPBH - Number of Pixels Per Block Horizontal
if obj.sicdmeta.ImageData.NumRows>8192
    NPPBV = 0; % (zero means "use NROWS")
else
    NPPBV = obj.sicdmeta.ImageData.NumRows;
end
fwriten(obj.FID, uint64(NPPBV), 4);  % NPPBV - Number of Pixels Per Block Vertical
fwriten(obj.FID, uint64(element_size),2);  % NBPP
fwriten(obj.FID, uint64(1), 3);  % IDLVL
fwriten(obj.FID, uint64(0), 3);  % IALVL
fwriten(obj.FID, uint64(0), 10); % ILOC - Image location offsets
fwriten(obj.FID, '1.0 ', 4);     % IMAG - image magnification
fwriten(obj.FID, uint64(0), 5);  % UDIDL
fwriten(obj.FID, uint64(0), 5);  % IXSHDL

% We're done so compute and return the number of bytes written
lengthWritten = ftell(obj.FID) - fileStartPosition;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////