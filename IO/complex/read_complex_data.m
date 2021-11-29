function [ complex_data, metadata ] = read_complex_data( filename, varargin )
%READ_COMPLEX_DATA Read complex SAR data from any format which is defined
%in the toolbox.
%
% [COMPLEX_DATA, METADATA] = READ_COMPLEX_DATA( FILENAME, DIMENSION_1_RANGE, ...
%    DIMENSION_2_RANGE, SUBSAMPLE, DECIMATION_FUNCTION, IMAGE_NUMBER)
%
% INPUTS:
%     FILENAME: Path and filename of complex data file.
%     DIMENSION_1_RANGE: 1x2 array [first element, last element].  If
%        reading a subimage, this is the range of data to read in the first
%        (azimuth) dimension.  Default is entire image range.  Can also use
%        an empty array ([]) to specify entire range.
%     DIMENSION_2_RANGE: 1x2 array [first element, last element].  If
%        reading a subimage, this is the range of data to read in the
%        second (range) dimension.  Default is entire image range.  Can
%        also use an empty array ([]) to specify entire range.
%     SUBSAMPLE:     1x2 array [subsample rate in first dimension,
%        subsample rate in second dimension].  Default is [1 1] (no
%        subsampling).
%     DECIMATION_FUNCTION: Decimation function, used if SUBSAMPLE is
%        defined to be greater than 1.  Options include 'none' (default)',
%        'max', 'mean', or any other function defined in Matlab.
%     IMAGE_NUMBER: For a file that contains multiple images, this is the
%        index of which file is to be read.  Default = 1.
% OUTPUTS:
%     COMPLEX_DATA: Array that contains the complex data from the AOI
%        selected.
%     METADATA: Structure that holds the metadata associated with the
%        file.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


if nargin < 6
    image_number=1;
else
    image_number=varargin{5};
end

file_reader_object=open_reader(filename);
if iscell(file_reader_object) % Multi-image file
    try % File must still be closed, even if read fails
        complex_data=file_reader_object{image_number}.read_chip(varargin{:});
        metadata=file_reader_object{image_number}.get_meta();
    catch
        for i=1:length(file_reader_object)
            file_reader_object{i}.close();
        end
        rethrow(lasterror);
    end
    for i=1:length(file_reader_object)
        file_reader_object{i}.close();
    end
else
    try % File must still be closed, even if read fails
        complex_data=file_reader_object.read_chip(varargin{:});
        metadata=file_reader_object.get_meta();
    catch
        file_reader_object.close();
        rethrow(lasterror);
    end
    file_reader_object.close();
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
