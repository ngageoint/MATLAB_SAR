function [ dataout ] = read_bip( fid, datasize, offset, ...
    datatype, bands, show_waitbar, varargin )
%READ_BIP Generic function for reading data band interleaved by pixel.
%
%   dataout = read_bip( fid, datasize, offset, datatype,...
%      bands, dim1range, dim2range, subsample )
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
%      BANDS:         Number of bands in data.  Default is 1.
%      SHOW_WAITBAR:  Determines whether read progress is displayed or not.
%         Default is true.
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
% Note: We write our own code, rather than use the MATLAB multibandread
% function, because the MATLAB multibandread function is unnacceptably
% SLOW!
%
% Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse function arguments and set defaults for undefined parameters
[dim1range,dim2range,subsample]=check_chipper_args(datasize, varargin{:});
% Sometimes datasize is passed in as a uint32.  If we don't cast to double,
% large files (>4 Gig) can result in overflows for some computed values.
datasize=double(datasize);
if((nargin<6)||isempty(show_waitbar)), show_waitbar=true; end;
if((nargin<5)||isempty(bands)), bands=1; end;
if((nargin<4)||isempty(datatype)), datatype='float32'; end;
[elementsize,memmapstr]=mdatatypeprops(datatype);
elementsize=elementsize*bands;
if((nargin<3)||isempty(offset)), offset=0; end;
offset=double(offset);
if length(offset)<2, offset(2)=0; end;
if((nargin<1)||(ftell(fid)<0))
    error('READ_BIP:InvalidArg','Invalid file identifier.');
end

%% Read data (region of interest only)
if fseek(fid,offset(1)+... % Beginning of data
    ((dim2range(1)-1)*(datasize(1)*elementsize+offset(2)))+... % Skip to first row of interest
    ((dim1range(1)-1)*elementsize),'bof') % Skip to first column of interest
    error('READ_BIP:EOF',['Attempt to read past end of file.  '...
     'File possibly corrupt.']);
end
dim1size=diff(dim1range)+1;
dim2size=diff(dim2range)+1;
if any(subsample>1) % Decimation required
    lengthdim2range=length(dim2range(1):subsample(2):dim2range(2));
    if verLessThan('matlab', '7.0')
        dataout=cast(zeros(length(dim1range(1):subsample(1):dim1range(2)),...
            lengthdim2range,bands),memmapstr);
    else
        dataout=zeros(length(dim1range(1):subsample(1):dim1range(2)),...
            lengthdim2range,bands,memmapstr);
    end
    if show_waitbar, h=waitbar(0,'Reading data...'); end
    for i=1:lengthdim2range
        single_line=fread(fid,double(bands*dim1size),[datatype '=>' memmapstr]);
        for j=1:bands % Pixel intervleaved
            dataout(:,i,j)=single_line(j:subsample(1)*bands:end);
        end
        fseek(fid,(((datasize(1)*elementsize)+offset(2))*(subsample(2)-1))+... % Skip unread rows
           ((datasize(1)-dim1size)*elementsize)+offset(2),'cof'); % Skip to beginning of dim2range
        if show_waitbar, waitbar(i/lengthdim2range,h); end;
    end
    if show_waitbar, close(h); end;
else % No decimation
    dataout=fread(fid,double([bands*dim1size,dim2size]),...
        [num2str(bands*dim1size) '*' datatype '=>' memmapstr],...
        double(offset(2)+(datasize(1)-dim1size)*elementsize));
    % Sometimes Matlab's FREAD function has an error
    % that results in an erronuous EOF flag when a skip value > 8192
    % (2^13) is specified.
    % if(feof(fid))  % check to see if the end of the file was overun
    %     error('READ_BIP:EOF',['Attempt to read past end of file.  '...
    %         'File possibly corrupt.']);
    % end
    if bands>1
        dataout=permute(reshape(dataout,[bands,dim1size,dim2size]),[2 3 1]); % Pixel intervleaved
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////