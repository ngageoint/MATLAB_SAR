function [ readerobj ] = open_cos_reader_noxml( filename, symmetry )
%OPEN_COS_READER_NOXML Intiates a reader object for TerraSAR-X COSAR file format.
%
% Does not look for the TerraSAR XML metadata file.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin<2
    symmetry=[0 0 1]; % COSAR written in range lines
end

%% Setup reader type
meta.native.cos=read_cos_meta(filename);
meta.ImageData.NumCols=uint32(meta.native.cos.az);
meta.ImageData.NumRows=uint32(meta.native.cos.rs);
datasize=[meta.native.cos.rs meta.native.cos.az];
datatype='int16';
complextype=true;
data_offset=[((datasize(1)+2)*4*4)+(2*4)... % Byte offset to first data sample 
             2*4]; % Byte spacing between rows
endian='b';
bands=1;

%% Open reader
readerobj=open_generic_reader(filename, datasize,...
    datatype, complextype, data_offset, endian, symmetry, bands);
readerobj.get_meta=@() meta;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////