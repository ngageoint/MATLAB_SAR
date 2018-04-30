function [ dataout ] = read_bib( fid, datasize, offset, ...
    datatype, bands_ip, blocksize, varargin )
%READ_BIP Generic function for reading data band interleaved by block.
%
%   dataout = read_bip( fid, datasize, offset, datatype,...
%      bands_ip, blocksize, dim1range, dim2range, subsample )
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
%      BANDS_IP:      Number of bands pixel interleaved within block.
%         Default is 1.
%      BLOCKSIZE:     1x2 array [size of block in first dimension, size of
%         block in second dimension].  Default is full image size.
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
% Warning: This code is extremely slow when reading large datasets, broken
% into many tiny blocks.
%
% Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse function arguments and set defaults for undefined parameters
[dim1range,dim2range,subsample]=check_chipper_args(datasize, varargin{:});
if((nargin<6)||isempty(blocksize)), blocksize=datasize; end;
if((nargin<5)||isempty(bands_ip)), bands_ip=1; end;
if((nargin<4)||isempty(datatype)), datatype='float32'; end;
[elementsize,memmapstr]=mdatatypeprops(datatype);
elementsize=elementsize*bands_ip;
if((nargin<3)||isempty(offset)), offset=0; end;
if((nargin<1)||(ftell(fid)<0))
    error('READ_BIP:InvalidArg','Invalid file identifier.');
end

%% Computer parameters necessary for quickly indexing into blocks
numTiles=ceil(datasize./blocksize);
tileBoundaries1=blocksize(1)*(1:numTiles(1))'; % Last index in the full image of each block (first dimension)
tileBoundaries1(:,2)=tileBoundaries1-blocksize(1)+1; % First index
tileBoundaries1(end,1)=datasize(1); % Might not be integer number of full blocks
tileBoundaries2=blocksize(2)*(1:numTiles(2))'; % Last index in the full image of each block (second dimension)
tileBoundaries2(:,2)=tileBoundaries2-blocksize(2)+1; % First index
tileBoundaries2(end,1)=datasize(2); % Might not be integer number of full blocks

regionIsInTile={dim1range(1)<=tileBoundaries1(:,1)&dim1range(2)>=tileBoundaries1(:,2),... % Dimension 1
    dim2range(1)<=tileBoundaries2(:,1)&dim2range(2)>=tileBoundaries2(:,2)}; % Dimension 2
tileRange={find(regionIsInTile{1},1,'first'):find(regionIsInTile{1},1,'last'),... % Dimension 1
    find(regionIsInTile{2},1,'first'):find(regionIsInTile{2},1,'last')}; % Dimension 2

tileSizes=ones(numTiles)*prod(blocksize)*elementsize; % First compute block sizes (in bytes)
tileOffsets=reshape(cumsum(tileSizes(:)),numTiles)-tileSizes+offset; % Array which shows offset in file to each block

% Setup output array
if verLessThan('matlab', '7.0')
    dataout=cast(zeros(length(dim1range(1):subsample(1):dim1range(2)),...
        length(dim2range(1):subsample(2):dim2range(2)),bands_ip),memmapstr);
else
    dataout=zeros(length(dim1range(1):subsample(1):dim1range(2)),...
        length(dim2range(1):subsample(2):dim2range(2)),bands_ip,memmapstr);
end

%% Read data (region of interest only)
use_waitbar=any(subsample>1)&&(length(tileRange{1})*length(tileRange{2}))>1;
if use_waitbar, h=waitbar(0,'Reading data...'); end;
xprogress=0;
for x=tileRange{1} % Iterate through tiles in first dimension
    yprogress=0;
    for y=tileRange{2} % Iterates through tiles in second dimension
        % Calculate indices relative to the full images of the region of
        % interest within only the current block
        roiIndicesFullImage={[max(dim1range(1), tileBoundaries1(x,2)),...
            min(dim1range(2), tileBoundaries1(x,1))],...
            [max(dim2range(1), tileBoundaries2(y,2)),...
            min(dim2range(2), tileBoundaries2(y,1))]};
        % Adjust to correct starting pixel within block taking subsampling into account
        roiIndicesFullImage{1}(1)=subsample(1)*...
            ceil((roiIndicesFullImage{1}(1)-dim1range(1))/subsample(1))+dim1range(1);
        roiIndicesFullImage{2}(1)=subsample(2)*...
            ceil((roiIndicesFullImage{2}(1)-dim2range(1))/subsample(2))+dim2range(1);
        % It is possible with severe subsampling to have blocks that have
        % no pixels to be read.
        if (diff(roiIndicesFullImage{1})>=0)&&(diff(roiIndicesFullImage{2})>=0)
            % Calculate indices within block
            roiIndicesTile={roiIndicesFullImage{1}-tileBoundaries1(x,2)+1,...
                roiIndicesFullImage{2}-tileBoundaries2(y,2)+1};
            % Read data from block
            block_data=read_bip( fid, blocksize, tileOffsets(x,y),...
                datatype, bands_ip, false,...
                roiIndicesTile{1}, roiIndicesTile{2}, subsample);
            if use_waitbar
                waitbar((xprogress*length(tileRange{2})+yprogress)/...
                    (length(tileRange{1})*length(tileRange{2})), h);
            end
            % Calculate indices within output images
            outputIndices={(roiIndicesFullImage{1}(1)-dim1range(1))/subsample(1),...
                (roiIndicesFullImage{2}(1)-dim2range(1))/subsample(2)};
            outputIndices{1}=outputIndices{1}+(1:length(roiIndicesTile{1}(1):subsample(1):roiIndicesTile{1}(2)));
            outputIndices{2}=outputIndices{2}+(1:length(roiIndicesTile{2}(1):subsample(2):roiIndicesTile{2}(2)));
            % Save to image
            dataout(outputIndices{:},:)=block_data;
        end
        yprogress=yprogress+1;
    end
    xprogress=xprogress+1;
end
if use_waitbar, close(h); end;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////