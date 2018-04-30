function [ mask ] = sicd_mask( sicd_meta, dim1_range, dim2_range )
%SICD_MASK Create a mask of valid data from SICD metadata
%
% SICD allows a ValidData construct that defines which pixels are valid
% values with a list of vertices of an enclosing polygon.  This functions
% converts that structure into a mask over any portion of the SICD image
% space.
%
% INPUTS:
%   sicd_meta  - required: SICD metadata structure.
%   dim1_range - optional: 1x2 array [first element, last element].  If
%        producing a mask for a subimage, this is the range of image area
%        in the first dimension (what is called "columns" in SICD-- but
%        not necessarily stored in that orientation in the MATLAB reader).
%        Default is entire image range. Can also use an empty array ([]) to
%        specify entire range.
%   dim2_range - optional: 1x2 array [first element, last element].  If
%        producing a mask for a subimage, this is the range of image area
%        in the second (range, called "rows" in SICD) dimension.  Default
%        is entire image range. Can also use an empty array ([]) to specify
%        entire range.
%
% OUTPUTS:
%   mask - Logical array over the area requested
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Default input parameters
if ~exist('row_range','var') || isempty(dim2_range)
    dim2_range = [1 sicd_meta.ImageData.FullImage.NumRows];
end
if ~exist('col_range','var') || isempty(dim1_range)
    dim1_range = [1 sicd_meta.ImageData.FullImage.NumCols];
end

if isfield(sicd_meta.ImageData,'ValidData') && ...
        isfield(sicd_meta.ImageData.ValidData,'Vertex')
    % Poly2mask has some peculiarities with including vertices (which are
    % document in the MATLAB help), so its possible that not all verticies
    % will be included in the resulting mask, but the boundary of the mask
    % will never be more than one pixel off of where it should be.
    mask = poly2mask(...
        [sicd_meta.ImageData.ValidData.Vertex.Row] + 1 ... % ValidData vertices are zero-based
           - dim2_range(1) + 1, ... % Relative to area of interest
        [sicd_meta.ImageData.ValidData.Vertex.Col] + 1 ... % ValidData vertices are zero-based
           - dim1_range(1) + 1, ... % Relative to area of interest
        diff(dim1_range)+1, diff(dim2_range)+1);
else % If valid data is not defined, assume all data is valid.
    mask = true(diff(dim1_range)+1, diff(dim2_range)+1);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////