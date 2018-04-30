function [ output_data, new_DeltaKCOAPoly ] = deskewmem( input_data, DeltaKCOAPoly, dim1_coords_m, dim2_coords_m, dim, fft_sgn )
%DESKEWMEM Performs deskew (centering of the spectrum on zero frequency) on a complex dataset
%    output_data = deskewmem(input_data, DeltaKCOAPoly, dim1_coords_m, dim2_coords_m, dim, fft_sgn)
%
%       Parameter name    Description
% 
%       input_data        Array of complex values to deskew
%       DeltaKCOAPoly     Polynomial that describes center of frequency
%                            support of data in the deskew dimension.
%                            Described in the SICD design and exploitation
%                            document.
%       dim1_coords_m     Coordinate of each "row" in dimension 1
%       dim2_coords_m     Coordinate of each "column" in dimension 2
%       dim               Dimension over which to perform deskew
%       fft_sgn           FFT sign required to transform data to spatial frequency domain
%       output_data       INPUT_DATA with proper phase adjustments to
%                            recenter frequency support at every point in
%                            the data
%       new_DeltaKCOAPoly Frequency support shift in the non-deskew
%                            dimension caused by the deskew.
%
% Implemented by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('fft_sgn','var')
    fft_sgn = -1; % Default, since it is more common
end

if dim==2, DeltaKCOAPoly=DeltaKCOAPoly.'; end;
% Integrate DeltaKCOA polynomial (in meters) to form new polynomial DeltaKCOAPoly_int
DeltaKCOAPoly_int=zeros(size(DeltaKCOAPoly)+[0 1]); % Integral has an extra row/column of zeros (constant term)
DeltaKCOAPoly_int(:,2:end)=DeltaKCOAPoly .* ...
    repmat(1./(1:size(DeltaKCOAPoly,2)),size(DeltaKCOAPoly,1),1); % Integrate
% New DeltaKCOAPoly in other dimension will be negative of the derivative
% of DeltaKCOAPoly_int in other dimension (assuming it was zero before).
new_DeltaKCOAPoly = -DeltaKCOAPoly_int(2:end,:) .* ...
    repmat((1:(size(DeltaKCOAPoly_int,1)-1)).',1,size(DeltaKCOAPoly_int,2));
if dim==2
    DeltaKCOAPoly_int=DeltaKCOAPoly_int.';
    new_DeltaKCOAPoly=new_DeltaKCOAPoly.';
end

% Apply phase adjustment from polynomial
output_data=input_data .* exp(1j * double(fft_sgn) * 2 * pi * ...
    sicd_polyval2d(DeltaKCOAPoly_int,dim1_coords_m,dim2_coords_m));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////