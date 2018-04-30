function [ readerobj ] = subset_reader(ro, dim1_limits, dim2_limits)
%SUBSET_READER Creates a reader with visibility into only part of a file
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Define object methods
readerobj = ro; % Copy most methods
readerobj.read_chip = @chip_fun; % Subset chipper
% Subset metadata
meta = readerobj.get_meta();
meta.ImageData.NumCols = diff(dim1_limits) + 1;
meta.ImageData.NumRows = diff(dim2_limits) + 1;
if isfield(meta.ImageData,'FirstCol')
    meta.ImageData.FirstCol = meta.ImageData.FirstCol + dim1_limits(1) - 1;
end
if isfield(meta.ImageData,'FirstRow')
    meta.ImageData.FirstRow = meta.ImageData.FirstRow + dim2_limits(1) - 1;
end
readerobj.get_meta = @() meta;

    % Check to see if indices are out of bounds of subset area and index
    % based on subset, not full image.
    function out = chip_fun(varargin)
        [ dim1range, dim2range, subsample ] =...
            check_chipper_args( [diff(dim1_limits) diff(dim2_limits)] + 1, varargin{:});
        out = ro.read_chip(dim1range + double(dim1_limits(1)) - 1, ...
            dim2range + double(dim2_limits(1)) - 1, subsample, varargin{4:end});
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////