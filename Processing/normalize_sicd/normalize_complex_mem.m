function [ output_data ] = normalize_complex_mem( input_data, ...
    DeltaKCOAPoly, dim1_coords_m, dim2_coords_m, ...
    weight_fun, oversample_rate, ...
    fftsign, dim )
%NORMALIZE_COMPLEX_MEM Normalize complex SAR data to bring to a consistent
%point in the processing chain
%
%    output_data = normalize_complex_mem(input_data, ...
%       DeltaKCOAPoly, dim1_coords_m, dim2_coords_m, ... % Deskew parameters
%       weight_fun, oversample_rate, ... % Deweighting paramters
%       fftsign, dim)
%
%       Parameter name    Description
% 
%       input_data        Array of complex values to deskew
%       DeltaKCOAPoly     Polynomial that describes center of frequency
%                            support of data in the processing dimension.
%                            Described in the SICD design and exploitation
%                            document.  DeltaKCOAPoly of 0 means no deskew.
%       dim1_coords_m     Coordinate of each "row" in dimension 1.  Used
%                            only for deskew.
%       dim2_coords_m     Coordinate of each "column" in dimension 2.  Used
%                            only for deskew.
%       weight_fun        Description of weighting applied.  Either a
%                            function handle to a function that takes a
%                            single argument (number of elements) and
%                            produces the weighting to apply, or a vector
%                            that is the weighting function sampled.  An
%                            empty weight_fun means to apply no
%                            deweighting.
%       oversample_rate   Amount of sampling beyond the ImpRespBW in the
%                            processing dimension.  (Nyquist sampling = 1).
%                            Used only for deweighting.
%       fftsign           FFT sign in the processing dimension (-1 or +1).
%       dim               Dimension over which to perform normalization
%       output_data       INPUT_DATA with normalization applied
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

output_data = input_data;
if any(DeltaKCOAPoly(:)~=0) % Do we need to center frequency support?
    output_data = deskewmem(input_data, DeltaKCOAPoly, dim1_coords_m, dim2_coords_m, dim, fftsign);
end
if fftsign==1 % -1 is the norm
    output_data = conj(output_data);
    % This transform will maintain the magnitude of the data and make the
    % IFFT behave like an FFT (and vice versa) in magnitude only:
    % abs(x)==abs(conj(x)) and
    % abs(fft(x))==abs(ifft(conj(x))), since fft(x)==conj(ifft(conj(x)))
end
if ~isempty(weight_fun) % Empty array means uniform weighting
    output_data = deweightmem(output_data, weight_fun, oversample_rate, dim);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////