function [complexdata, userdata] = read_sio(filename, varargin)
%READ_SIO Read SIO formatted complex data
%   complexdata = read_sio( filename, rowrange, colrange, subsample, decimationfun )
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
%      USERDATA: SIO metadata passed as "user data" in SIO header.
%         (Typically this is an empty structure.)
%
%   Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%read_only, big-endian, all non-header data 32-bit float
[ ihdr, endian, data_offset, userdata ] = read_sio_meta( filename );
[datatype,complexbool,freadtype]=siotype2matlab(ihdr(4),ihdr(5));

fid = fopen(filename,'r',endian); 
complexdata=read_complex(fid,[ihdr(3) ihdr(2)],data_offset,freadtype,complexbool,varargin{:});
fclose(fid);

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////