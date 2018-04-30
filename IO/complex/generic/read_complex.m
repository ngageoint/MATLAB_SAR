function [ complexdata ] = read_complex( fid, datasize, offset, ...
    datatype, complexbool, varargin )
%READ_COMPLEX Generic function for reading binary complex data.
%
%   Assumes i and q components interleaved (adjacent for a single pixel).
%   Uses regular fread functions to read in data.
%
%   complexdata = read_complex( fid, datasize, offset, datatype,...
%      complexbool, dim1range, dim2range, subsample, decimationfun )
%
%   INPUTS:
%      FID:           MATLAB file identifier.  Must refer to a file that is
%         open for reading.
%      DATASIZE:      1x2 array [number of elements in first dimension,
%         number of elements in the second dimension].  Total number of
%         elements in each dimension of the data, with dimension ordered as
%         there are written in the file.
%      OFFSET:        Index (in bytes) from the beginning of the file to
%         the beginning of the data.  Default is 0 (beginning of file).
%      DATATYPE:      MATLAB string for specifying binary data precision
%         ('uint8','float32',etc.) See FREAD documentation for full list.
%         Default is 'float32'.
%      COMPLEXBOOL:   Boolean whether data is complex (true) or real
%      DIM1RANGE:     1x2 array [first element, last element].  If reading
%         a subimage, this is the range of data to read, in the first
%         dimension (as written in the file).  Default is entire image
%         range.  Can also use an empty array ([]) to specify entire range.
%      DIM2RANGE:     1x2 array [first element, last element].   Default is
%         entire image range.  Can also use an empty array ([]) to specify
%         entire range.
%      SUBSAMPLE:     1x2 array [subsample rate in first dimension,
%         subsample rate in second dimension].  Default is [1 1] (no
%         subsampling).
%
%   OUTPUTS:
%      COMPLEXDATA:   Array of complex data values, of data type single
%         complex.  Is oriented so that the first dimension in the file
%         order matches the first dimension in MATLAB (vertical when
%         display in imshow, imagesc, etc.), so this may need to be
%         transposed for display for many image types.
%
% Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if complexbool
    complexdata=read_bip(fid, datasize, offset, datatype, 2, true, varargin{:});
    complexdata=complex(complexdata(:,:,1:2:end),complexdata(:,:,2:2:end));
else
    complexdata=read_bip(fid, datasize, offset, datatype, 1, true, varargin{:});
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////