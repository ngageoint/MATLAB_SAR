function  [ output_image, meta ]= read_gff( filename, varargin )
%READ_GFF Sandia GSAT Image File Format
%   complexdata = read_gff( filename, rowrange, colrange, subsample, decimationfun )
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
%      META:          Structure containing metadata
%
%   Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Read metadata
native_meta = read_gff_meta(filename);
meta=meta2sicd_gff(native_meta);
meta.native.gff=native_meta;
if native_meta.RowMajor
   datasize = [native_meta.RgCnt native_meta.AzCnt];
else
   datasize = [native_meta.AzCnt native_meta.RgCnt];
end
switch native_meta.ImageType
    case 0
        datatype='uchar';
    case 1
        if native_meta.BytesPerPixel == 4 % 2 bytes phase, 2 bytes magnitude
            datatype='uint16';
        else % 8 bytes total, 4 bytes phase, 2 bytes magnitude
            datatype='uint32';
        end
    case 2
        datatype='float32';
end
if native_meta.Endian
    endian='b';
else
    endian='l';
end

% Parse input parameters; reverse indices since image 180 degrees rotated
% from shadows down orientations
newvarargin=varargin;
if nargin>1
    rowrange=datasize(1)-varargin{1}(end:-1:1)+1;
    if native_meta.RowMajor
        newvarargin{2}=rowrange;
        newvarargin{1}=[];
    else
        newvarargin{1}=rowrange;
    end
end
if nargin>2
    colrange=datasize(2)-varargin{2}(end:-1:1)+1;
    if native_meta.RowMajor
        newvarargin{1}=colrange;
    else
        newvarargin{2}=colrange;
    end
end
if nargin>3
    if native_meta.RowMajor
        newvarargin{3}=varargin{3}(end:-1:1);
    end
end

% Read complex data
fid = fopen(filename,'r',endian);
output_image = read_complex(fid,datasize,native_meta.Length,...
    datatype,logical(native_meta.ImageType),newvarargin{:});
fclose(fid);

if native_meta.ImageType==1 % Int types are phase-magnitude, rather than real-imag, which read_complex expects
    output_image=double(imag(output_image)).*exp(1i * double(real(output_image)) * 2*pi/(2^16));
end
if native_meta.RowMajor
    output_image=output_image.';
end
output_image=rot90(output_image,2); % GFF stored shadows up; switch to shadows down

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////