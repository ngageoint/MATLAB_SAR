function [lats, lons] = project_image(beginPoint, imageSize, skip, meta, varargin)
%PROJECT_IMAGE "Projects" the coordinates from each pixel in the supplied
% slant-plane image to the ground-plane.
%
% Inputs:
%       beginPoint:       The starting point (row,column) of chip to project
%       imageSize:        The size of the image chip to project (rows, cols)
%       skip:             The "stride" of pixels to extract within the
%                         specified region.  A skip of 1 will take every
%                         pixel, a skip of 2 will take every other pixel,
%                         etc.
%       meta:             Metadata for the image (as loaded by the 
%                         open_reader "get_meta" function)
%
% Outputs:
%       lats:             The latitudes locations for each image pixel
%       lons:             The longitude locations for each image pixel
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

    [points_col, points_row] = meshgrid(beginPoint(1) + (0:skip(1):(imageSize(1)-1)), ...
                                       beginPoint(2) + (0:skip(2):(imageSize(2)-1)));
    pos_lla = point_slant_to_ground([points_row(:)'; points_col(:)'], meta, varargin{:});
    outputImageSize=size(points_col);
    lats = reshape(pos_lla(1,:),outputImageSize(1),outputImageSize(2))';
    lons = reshape(pos_lla(2,:),outputImageSize(1),outputImageSize(2))';

end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
