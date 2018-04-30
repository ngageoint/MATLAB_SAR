function [ readerobj ] = open_ctiff_reader_noxml( filename, symmetry )
%OPEN_CTIFF_READER Read a TIFF that uses a complex SampleFormat
%
% Unfortunately, we have to write our own code here since the TIFF library
% that MATLAB uses does not support a complex SampleFormat.  Actually, this
% reader works for other TIFFs as well, but for those we can use the
% builtin MATLAB functions.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin<2
    symmetry=[0 0 0];
end

%% Get TIFF metadata
tags = read_tiff_tags(filename);

%% Error check
% This function is intentionally only limited capability, purely to
% overcome the fact that the MATLAB TIFF reader doesn't handle complex data
% types.  None of TIFF files with complex data that we have ever seen have
% fallen into one of these cases checked for below.  If we find one, we may
% have to extend the reader to handle that case.  These checks will at
% least warn a user to explain why the reader is not working.
fixed_rowlength = all(tags{1}.StripByteCounts==tags{1}.StripByteCounts(1));
contiguous = all((tags{1}.StripOffsets(1:(end-1)) + ...
    tags{1}.StripByteCounts(1:(end-1))) == tags{1}.StripOffsets(2:end));
if ~fixed_rowlength || ~contiguous
    error('OPEN_CTIFF_READER_NOXML:NONCONTIGUOUS_DATA','This TIFF reader only handles contiguous data.');
end
if isfield(tags{1},'Compression') && tags{1}.Compression~=1
    error('OPEN_CTIFF_READER_NOXML:COMPRESSED_DATA','This TIFF reader does not handle compressed data.');
end
if isfield(tags{1},'FillOrder') && tags{1}.FillOrder~=1
    error('OPEN_CTIFF_READER_NOXML:COMPRESSED_DATA','This TIFF reader does not handle bit-reversed data.');
end
% Reader currently only handle pixel-interleaved ("chunky") data
if isfield(tags{1},'PlanarConfiguration') && tags{1}.PlanarConfiguration~=1
    error('OPEN_CTIFF_READER_NOXML:COMPRESSED_DATA','This TIFF reader does not handle "planar" format (band-interleaved) data.');
end
if isfield(tags{1},'ExtraSamples') && ~isempty(tags{1}.ExtraSamples)
    error('OPEN_CTIFF_READER_NOXML:EXTRA_SAMPLES','This TIFF reader does not handle extra samples per pixel.');
end
if isfield(tags{1},'TileOffsets')
    error('OPEN_CTIFF_READER_NOXML:COMPRESSED_DATA','This TIFF reader does not handle tiled data.');
end

%% Testing pixel reading here
% For now, we just use first image.  Later if we run into TIFF files with
% more than one image, we will have to include all of them and return a
% cell array of reader objects.
image_n = 1;
% Setup SICD metadata stub
meta.native.tiff=tags{image_n};
if symmetry(3)
    meta.ImageData.NumCols=uint32(tags{image_n}.ImageLength);
    meta.ImageData.NumRows=uint32(tags{image_n}.ImageWidth);
else
    meta.ImageData.NumCols=uint32(tags{image_n}.ImageWidth);
    meta.ImageData.NumRows=uint32(tags{image_n}.ImageLength);
end
datasize=double([tags{image_n}.ImageWidth tags{image_n}.ImageLength]);
% Compute parameters necessary for reading pixel data
% Determine TIFF endianness
fid = fopen(filename,'r');
endian = fread(fid,2,'uint8=>char').';
fclose(fid);
switch endian
    case 'II'
        endian = 'l';
    case 'MM'
        endian = 'b';
    % No need for otherwise since this would have already been checked in
    % the read_tiff_tags() call above.
end
if isfield(tags{image_n},'SampleFormat')
    switch tags{image_n}.SampleFormat(1) % Assumes same format for all bands
        case 1
            datatype = 'uint';
            complextype = 0;
        case 2
            datatype = 'int';
            complextype = 0;
        case 3
            datatype = 'float';
            complextype = 0;
        case 5
            datatype = 'int';
            complextype = 1;
        case 6
            datatype = 'float';
            complextype = 1;
        otherwise
            error('Unrecognized datatype.');
    end
else
    datatype = 'uint';
end
bits_per_sample = tags{image_n}.BitsPerSample(1) / (complextype+1); % Assumes same format for all bands
datatype = [datatype num2str(bits_per_sample)];
bands = double(tags{image_n}.SamplesPerPixel);
% Initialize reader object
readerobj=open_generic_reader(filename, datasize,...
    datatype, complextype, double(tags{image_n}.StripOffsets(1)), endian, symmetry, bands);
readerobj.get_meta=@() meta;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////