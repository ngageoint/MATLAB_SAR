function [ output ] = sicd_polyval2d( poly_coefs, dim1_ind, dim2_ind, sicd_meta )
%SICD_POLYVAL2D Evaluate a SICD field of type 2D_POLY
%
% MATLAB doesn't have a builtin function for evaluating 2D polynomials, so
% this is the function we use for evaluating the 2D polynomials in SICD
% over a grid of input values. This function also accepts image indices
% (which is more natural in MATLAB) rather than image coordinates (which is
% how SICD 2D_POLYs are generally defined).
%
% INPUTS:
%   poly_coefs - required: A SICD field of type 2D_POLY.
%   dim1_ind   - required: Vector of value(s) in the first dimension
%                (azimuth).
%   dim2_ind   - required: Vector of value(s) in the second dimension
%                (range).  Can be a vector.
%   sicd_meta  - optional: SICD metadata structure.  If this is passed,
%                dim1_ind and dim2_ind will be converted from image indices
%                to image coordinates (meters from SCP).  Otherwise values
%                are used directly with no conversion.
%
% OUTPUTS:
%   output - Array with the result values of the evaluated polynomial.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% SICD uses the first dimension of its polynomials to refer to the "row"
% index (range) and the second dimension to refer to the "column" index
% (azimuth).  In our MATLAB framework, we store data the opposite way; that
% is, the first index (rows in MATLAB) is used for azimuth and the second
% index (columns in MATLAB) is used for range.
poly_coefs = poly_coefs.'; % Swap rows/columns

if exist('sicd_meta','var') && ~isempty(sicd_meta)
    % Convert image indices (irow and icol in the notation used in the SICD
    % spec) to image coordinates (xrow and ycol in the notation used in the
    % SICD spec) in meters from SCP.
    dim1_vals = (double(dim1_ind-1) + double(sicd_meta.ImageData.FirstCol) - ...
        double(sicd_meta.ImageData.SCPPixel.Col)) * sicd_meta.Grid.Col.SS;
    dim2_vals = (double(dim2_ind-1) + double(sicd_meta.ImageData.FirstRow) - ...
        double(sicd_meta.ImageData.SCPPixel.Row)) * sicd_meta.Grid.Row.SS;
else % No conversion requested
    dim1_vals = dim1_ind;
    dim2_vals = dim2_ind;
end

% Evaluate polynomial
output=zeros(length(dim1_vals),length(dim2_vals));
[x_2d,y_2d]=ndgrid(dim1_vals,dim2_vals);
for i=1:size(poly_coefs,1)
    for j=1:size(poly_coefs,2)
        output=output+poly_coefs(i,j).*(x_2d.^(i-1)).*(y_2d.^(j-1));
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////