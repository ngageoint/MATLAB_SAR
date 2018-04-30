function [ output_image, info ] = read_tiff( filename, varargin )
%READ_TIFF Read in TIFF file
%   complexdata = read_tiff( filename, rowrange, colrange, subsample, decimationfun )
%
%   Works for complex RADARSAT-2 GeoTIFF data
%
%   INPUTS:
%      FILENAME:      Input file to read.
%      ROWRANGE:      1x2 array [first row, last row], row #1 is top.
%         Default is entire image range.  Can also use an empty array ([])
%         to specify entire row range.
%      COLRANGE:      1x2 array [first column, last column], column #1 is
%         left side.   Default is entire image range.  Can also use an
%         empty array ([]) to specify entire column range.
%      SUBSAMPLE:     1x2 array [subsample rate in first dimension,
%         subsample ratein second dimension].  Default is [1 1] (no
%         subsampling).
%      Note: There are no DECIMATION options as in most of the other READ_*
%      functions.
%
%   OUTPUTS:
%      COMPLEXDATA:   Array of complex data values of data type single
%         complex.
%
%   Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if(nargin<4)
    subsample=[1 1];
else
    subsample=varargin{3};
end
info=imfinfo(filename);
if(nargin<3)||isempty(varargin{2})
    colrange=[1 info.Width];
else
    colrange=varargin{2};
end
if(nargin<2)||isempty(varargin{1})
    rowrange=[1 info.Height];
else
    rowrange=varargin{1};
end

multiband_image = imread(filename,'tiff', 'PixelRegion', ...
    {[rowrange(1) subsample(1) rowrange(2)],...
     [colrange(1) subsample(2) colrange(2)]});
output_image=complex(multiband_image(:,:,1),multiband_image(:,:,2));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////