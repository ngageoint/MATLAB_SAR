function [ readerobj ] = open_generic_reader(filename, datasize,...
    datatype, complextype, data_offset, endian, symmetry, bands_ip,...
    blocksize)
%OPEN_GENERIC_READER Generic reader object for reading pixel interleaved
% data from a flat file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if (nargin<7||isempty(symmetry)), symmetry=[0 0 0]; end; % Assume store azimuth-major order, viewed from above
if (nargin<8||isempty(bands_ip)), bands_ip=1; end; % Assume single band
if (nargin<9||isempty(blocksize)), blocksize=datasize; end; % Assume no blocking

[chipper_function, close_function]=generic_chipper( filename, datasize,...
    datatype, complextype, data_offset, endian, symmetry, bands_ip, blocksize);
if symmetry(3), datasize=datasize(end:-1:1); end;
readerobj = chipfun2readerobj(chipper_function, datasize);
readerobj.close=close_function;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////