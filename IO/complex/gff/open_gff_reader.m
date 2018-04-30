function [ readerobj ]= open_gff_reader( filename )
% OPEN_GFF_READER Initiates a reader object for GFF file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse metadata
native_meta=read_gff_meta(filename);
meta=meta2sicd_gff(native_meta);
meta.native.gff=native_meta;
datasize=[meta.ImageData.NumCols meta.ImageData.NumRows];
if native_meta.RowMajor
    datasize=datasize(end:-1:1);
end
switch native_meta.ImageType
    case 0  % Magnitude only
        datatype='uchar';
        complextype=1;
    case 1 % Phase/magnitude
        if native_meta.BytesPerPixel == 4 % 2 bytes phase, 2 bytes magnitude
            datatype='uint16';
        else % 8 bytes total, 4 bytes phase, 4 bytes magnitude
            datatype='uint32';
        end
        complextype=@(x) double(x(:,:,2)).*exp(2 * pi * 1i * double(x(:,:,1)) /...
            double(intmax(datatype)+1));
    case 2 % I/Q float
        datatype='float32';
        complextype=2;
end
if native_meta.Endian
    endian='b';
else
    endian='l';
end
data_offset=native_meta.Length;
% GFF stored shadows up; switch to shadows down
if native_meta.RowMajor
    symmetry=[1 1 1];
else
    symmetry=[1 1 0];
end
bands=1;

%% Build object
readerobj=open_generic_reader(filename, datasize,...
    datatype, complextype, data_offset, endian, symmetry, bands);
readerobj.get_meta=@() meta;

end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
