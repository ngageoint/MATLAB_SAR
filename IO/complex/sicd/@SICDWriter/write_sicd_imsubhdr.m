function write_sicd_imsubhdr(obj, n)
%WRITE_SICD_IMSUBHDR Writes the NITF image segment sub-header for SICD file
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
% Written by: Tom Krauss and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Always writing version 2.1 NITF

fwriten(obj.FID, 'IM', 2);
if obj.NumIS==1
    IID1 = 'SICD000';
else
    IID1 = ['SICD' num2str(n,'%03d')];
end
fwriten(obj.FID, IID1, 10); % IID1

if isfield(obj.sicdmeta,'ImageCreation')&&isfield(obj.sicdmeta.ImageCreation,'DateTime')
    fdt = datestr(obj.sicdmeta.ImageCreation.DateTime,'yyyymmddHHMMSS');
else
    fdt = '';
end
fwriten(obj.FID, fdt, 14); % IDATIM

fwriten(obj.FID, '', 17); % TGTID

if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'CoreName')
    IID2 = obj.sicdmeta.CollectionInfo.CoreName;
else
    IID2 = 'Unknown';
end
fwriten(obj.FID, ['SICD: ' IID2], 80); % IID2

obj.write_sicd_security_tags();  % Security tags go here
fwriten(obj.FID, '0', 1);  % ENCRYPT
if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'CollectorName')
    ISORCE = obj.sicdmeta.CollectionInfo.CollectorName;
else
    ISORCE = '';
end
fwriten(obj.FID, ['SICD: ' ISORCE], 42); % ISORCE
fwriten(obj.FID, uint64(obj.NumRowsIS(n)), 8); % NROWS
fwriten(obj.FID, uint64(obj.sicdmeta.ImageData.NumCols), 8); % NCOLS

switch obj.sicdmeta.ImageData.PixelType
    case 'RE32F_IM32F'
        PVTYPE = 'R';
        ABPP = 32;
        ISUBCAT1 = 'I';
        ISUBCAT2 = 'Q';
    case 'RE16I_IM16I'
        PVTYPE = 'SI';
        ABPP = 16;
        ISUBCAT1 = 'I';
        ISUBCAT2 = 'Q';
    case 'AMP8I_PHS8I'
        PVTYPE = 'INT';
        ABPP = 8;
        ISUBCAT1 = 'M';
        ISUBCAT2 = 'P';
end
fwriten(obj.FID, PVTYPE, 3); % PVTYPE (here signed integer)
fwriten(obj.FID, 'NODISPLY',8); % IREP
fwriten(obj.FID, 'SAR', 8); % ICAT
fwriten(obj.FID, uint64(ABPP),2);  % ABPP
fwriten(obj.FID, 'R', 1);  % PJUST
fwriten(obj.FID, 'G', 1);  % ICORDS
if isfield(obj.sicdmeta,'GeoData')&&isfield(obj.sicdmeta.GeoData,'ImageCorners')
    geo_format = {'num_units',3,'include_symbols',false};
    % TODO: Adjust to handle multiple image segments
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
fwriten(obj.FID, geolo, 60); % IGEOLO

fwriten(obj.FID, uint64(0), 1);  % NICOM (# image comments)
fwriten(obj.FID, 'NC', 2);  % IC (No image compression)
fwriten(obj.FID, uint64(2), 1);  % NBANDS

% 'Band' 1: in phase
fwriten(obj.FID, '',  2); % IREPBAND1
fwriten(obj.FID, ISUBCAT1, 6); % ISUBCAT1
fwriten(obj.FID, 'N', 1); % IFC1
fwriten(obj.FID, '',  3); % IMFLT1
fwriten(obj.FID, '0', 1); % NLUTS1

% 'Band' 2: quadrature
fwriten(obj.FID, '',  2); % IREPBAND2
fwriten(obj.FID, ISUBCAT2, 6); % ISUBCAT2
fwriten(obj.FID, 'N', 1); % IFC2
fwriten(obj.FID, '',  3); % IMFLT2
fwriten(obj.FID, '0', 1); % NLUTS2

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
if obj.NumRowsIS(n)>8192
    NPPBV = 0; % (zero means "use NROWS")
else
    NPPBV = obj.NumRowsIS(n);
end
fwriten(obj.FID, uint64(NPPBV), 4);  % NPPBV - Number of Pixels Per Block Vertical
fwriten(obj.FID, uint64(ABPP),2);  % NBPP
fwriten(obj.FID, uint64(n), 3);  % IDLVL
fwriten(obj.FID, uint64(n-1), 3);  % IALVL
if n==1
    fwriten(obj.FID, uint64(0), 5); % ILOC (row) - Image location offsets
else
    fwriten(obj.FID, uint64(obj.NumRowsIS(n-1)), 5); % ILOC (row) - Image location offsets
end
fwriten(obj.FID, uint64(0), 5); % ILOC (column) - Image location offsets
fwriten(obj.FID, '1.0 ', 4);     % IMAG - image magnification
fwriten(obj.FID, uint64(0), 5);  % UDIDL - Generally not used for SICD
fwriten(obj.FID, uint64(0), 5);  % IXSHDL - Generally not used for SICD

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////