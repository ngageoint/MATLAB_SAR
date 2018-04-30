function [complexdata, meta] = read_cos(filename, rowrange, colrange, subsample, decimationfun)
% Read TerraSAR-X COSAR file format.
%
%   [ complexdata, metadata ] =...
%      read_cos( filename, rowrange, colrange, subsample, decimationfun )
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
%      DECIMATIONFUN: Decimation function, used if SUBSAMPLE is defined to
%         be greater than 1.  Options include 'none' (default)', 'max',
%         'mean', or any other function defined in Matlab.
%
%   OUTPUTS:
%      COMPLEXDATA:   Array of complex data values of data type single
%         complex.
%      METADATA:      Matlab struct that contains the metadata within the
%         file.
%  
%   DEPENDENCIES:
%      MATLAB Image Processing Toolbox: Is required if subsampling option
%      is used.  blkproc and/or im2col are the only functions used from
%      this toolbox.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Read metadata file
meta=read_cos_meta(filename);

% Open files to read data
fid = fopen(filename,'r','b'); 

% COSAR stores in column(range)-major order so switch order of arguments
% to what READ_COMPLEX expects (row-major)
if(nargin<5), decimationfun='none'; end;
if(nargin<4), subsample=[1 1]; end;
varargin={rowrange+2, colrange+1, subsample(end:-1:1), decimationfun};

% Read binary complex data
complexdata=read_complex(fid,[meta.rs meta.az]+[2 1],0,...
    'int16',true,varargin{:});

% Reorient data to be energy-from-top
complexdata = complexdata.';

% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////