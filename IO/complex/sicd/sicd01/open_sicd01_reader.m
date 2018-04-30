function [ readerobj ] = open_sicd01_reader( filename )
%OPEN_SICD01_READER Intiates a reader object for Sensor Independent Complex
%Data (SICD) file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup reader
meta.native.sicd01=read_sicd01_meta(filename);
meta.ImageData.NumCols=meta.native.sicd01.header.data_columns;
meta.ImageData.NumRows=meta.native.sicd01.header.data_rows;
datasize=[meta.ImageData.NumCols meta.ImageData.NumRows];
datatype='float32';
complexbool=true;
data_offset=meta.native.sicd01.header_length_in_bytes+meta.native.sicd01.header.xml_content_length+4;
endian='b';

%% Build object
readerobj=open_generic_reader(filename, datasize, datatype, complexbool,...
    data_offset, endian);
readerobj.get_meta=@() meta;
readerobj.get_thumb=@() read_sicd01_thumb(filename, meta.native.sicd01);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////