function [ readerobj ] = open_reader( filename )
%OPEN_READER Open a generic reader for multiple formats of complex imagery
%
% Returns a reader object structure with the following methods:
%
% readerobj.read_chip(dimension_1_range, dimension_2_range, subsample, decimation_function)
%      DIMENSION_1_RANGE:     1x2 array [first element, last element].  If
%         reading a subimage, this is the range of data to read, in the
%         first (azimuth) dimension.  Default is entire image range.  Can
%         also use an empty array ([]) to specify entire range.
%      DIMENSION_2_RANGE:     1x2 array [first element, last element].  If
%         reading a subimage, this is the range of data to read, in the
%         second (range) dimension. Default is entire image range.  Can
%         also use an empty array ([]) to specify entire range.
%      SUBSAMPLE:     1x2 array [subsample rate in first dimension,
%         subsample rate in second dimension].  Default is [1 1] (no
%         subsampling).
%      DECIMATION_FUNCTION: Decimation function, used if SUBSAMPLE is
%         defined to be greater than 1.  Options include 'none' (default)',
%          'max', 'mean', or any other function defined in Matlab.
% readerobj.get_meta() % Returns a MATLAB structure with all metadata
%      contained in the file
% readerobj.close() % Should always be run after all reading is finished
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


%% Determine format and open
format_string=guess_complex_format(filename);
if isempty(format_string)
    if ischar(filename) && ~exist(filename,'file')
        error('OPEN_READER:FileDoesNotExist',['File ''' filename ''' does not exist.']);
    else
        error('OPEN_READER:UnrecognizedFormat','Unable to determine complex image format.');
    end
end
readerobj = feval(['open_' format_string '_reader'], filename);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////