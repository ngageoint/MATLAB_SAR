function [ dataout ] = read_bip_mm( fileobj, endian, varargin )
%READ_BIP_MM Generic function for reading band-interleaved-by-pixel data using memory mapping
%
%   dataout = read_bip_mm( fileobj, endian, complexbool, dim1range, dim2range, subsample )
%
%   INPUTS:
%      FILEOBJ:       MATLAB memmapfile object with data in fileobj.data.val
%      ENDIAN:        Big-endian ('b') or little-endian ('l') data.
%         Default is 'b'.
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
%      dataout:   Array of complex data values, of data type single
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

%% Parse arguments
datasize=size(fileobj.data.val); datasize=datasize([2 3]); % Ignore number of bands
[dim1rng,dim2rng,subsamp]=check_chipper_args(datasize, varargin{:});
if((nargin<3)||isempty(endian)), endian='b'; end;
[comptype,varsize,compend]=computer;
swap=(upper(endian)~=upper(compend));

%% Chip file (in row and column, but read all bands)
dataout=fileobj.data.val(:,dim1rng(1):subsamp(1):dim1rng(2),dim2rng(1):subsamp(2):dim2rng(2));
dataout=permute(dataout,[2 3 1]); % Band interleaved by pixel.  Puts bands in 3rd dimension
if(swap), dataout=swapbytes(dataout); end % Is native byte order same as file byte order?

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////