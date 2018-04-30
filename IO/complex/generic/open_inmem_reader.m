function [ readerobj ] = open_inmem_reader( filename )
%OPEN_INMEM_READER Intiates a reader object for arrays of data in memory
%
% Typically the open_reader framework handles complex data files.  This
% function handles the special case when an array of data in memory is
% passed, rather than a filename.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if isnumeric(filename)
    data = filename;
    meta.ImageData.NumRows = size(data,2);
    meta.ImageData.NumCols = size(data,1);
else % Cell array with data and SICD metadata structure
    data = filename{1};
    meta = filename{2};
end

readerobj.read_chip = @chipper;
readerobj.get_meta = @(varargin) meta;
readerobj.close = @(varargin) deal(); % Do nothing

    function out = chipper(varargin)
        [dim1range,dim2range,subsample]=check_chipper_args(size(data), varargin{:});
        if nargin>=4 && any(subsample>1) && ...
                ~strncmpi(varargin{4},'none',4) && ~isempty(varargin{4})
            out = data(dim1range(1):dim1range(2),dim2range(1):dim2range(2)); % Before subsampling
            out = reshape(feval(varargin{4},im2col(out,subsample,'distinct')), ...
                ceil(size(out)./subsample)); % Apply decimation
        else
            out = data(dim1range(1):subsample(1):dim1range(2),...
                dim2range(1):subsample(2):dim2range(2));
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////