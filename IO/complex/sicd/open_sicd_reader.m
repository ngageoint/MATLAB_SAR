function [ readerobj ] = open_sicd_reader( filename )
%OPEN_SICD_READER Intiates a reader object for Sensor Independent Complex
%Data (SICD) file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup reader type
[SICD_meta, NITF_meta]=read_sicd_meta(filename);
if isfield(NITF_meta,'matlab_sar')
    SICD_meta.native.nitf = NITF_meta.matlab_sar;
end

%% NITF container header data
data_offset=NITF_meta.minimal.imageSegmentOffsets;
datasize=double([NITF_meta.minimal.imageSegmentColumns, NITF_meta.minimal.imageSegmentRows]);

%% SICD XML header data
switch SICD_meta.ImageData.PixelType
    case 'RE32F_IM32F'
        datatype='float32';
        complexbool=true;
    case 'RE16I_IM16I'
        datatype='int16';
        complexbool=true;
    case 'AMP8I_PHS8I'
        datatype='uint8';
        if ~isfield(SICD_meta.ImageData,'AmpTable')
            % If the optional table isn't there, value should be used directly.
            SICD_meta.ImageData.AmpTable = 0:255;
        end
        complexbool=@(x) SICD_meta.ImageData.AmpTable(x(:,:,1)+1).* ...
            exp(2*pi*1i*single(x(:,:,2))/256);
    otherwise
        error('OPEN_SICD_READER:INVALID_PIXEL_TYPE','Invalid pixel type.');
end
endian='b';
symmetry=[0 0 0];
bands=1;

%% Build object
readerobj=multisegment_reader(filename, datasize, datatype, complexbool,...
    data_offset, endian, symmetry, bands);
readerobj.get_meta=@() SICD_meta;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////