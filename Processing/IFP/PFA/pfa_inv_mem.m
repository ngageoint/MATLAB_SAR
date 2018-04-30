function [resampled_data, k_a, k_r] = pfa_inv_mem(complex_data, center_freq, ss, imp_resp_bw, fft_size, fft_sgn)
%PFA_INV_MEM Inverse polar format IFP algorithm on an array of complex data
%
% Inputs:
%     complex_data:   Complex SAR data.  Arranged range-down,
%                     view-from-above.  Although this is transposed from
%                     the typical MATLAB display convention (range is along
%                     the second dimension).
%     center_freq:    Center frequency of rectangular data in cycles/meter.
%                     This is equivalent to Grid.Row.KCtr in SICD.
%     ss:             Sample spacing across [col, row].
%     imp_resp_bw:    Optional.  Impulse response bandwidth across
%                     [col,row] as defined in SICD.  Default is no
%                     zeropad.
%     fft_size:       Optional.  Size of resampled data.  Default is
%                     COMPLEX_DATA size.
%     fft_sgn:        Optional.  Sign of the FFT to transform into the
%                     spatial frequency domain.  Default is -1.
%
% Outputs:
%     resampled_data: Data uniformly sampled in the k_a/k_r domain.
%     k_a:            Bounds of data in polar angle (radians).
%     k_r:            Bounds of data in frequency (cycles/meter).
%
% This function is the inverse of pfa_mem, useful for transforming and
% resampling data to a frequency/polar angle grid.
%
% Data is assumed to have its frequency support centered.
%
% Authors: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Default values for input parameters
if ~exist('imp_resp_bw','var')
    imp_resp_bw = 1./ss; % No zeropad
end
if ~exist('fft_size','var')
    fft_size = size(complex_data);
end
if ~exist('fft_sgn','var')
    fft_sgn = -1; % More common
end

%% Transform data from image domain to spatial frequency domain
if fft_sgn<0
    fft_fun = @fft2;
else
    fft_fun = @ifft2;
end
complex_data = fftshift(fft_fun(complex_data,fft_size(1),fft_size(2)));
complex_data = complex_data/sqrt(prod(fft_size));

%% Compute new coordinates onto which to interpolate
% We have to compute the bounds for the inscribed annulus that just fits in
% our rectangularly sampled input data without zeropad.
k_v = center_freq + (imp_resp_bw(2) * [-1 1] / 2); % Coordinates of nonzeropad area
k_u = imp_resp_bw(1) * [-1 1] / 2; % Coordinates of nonzeropad area
k_a = asin(k_u/k_v(2)); % angular extents of the inscribed annulus in radians
k_r = [k_v(1)./cos(k_a(1)), k_v(2)]; % radial extents of the inscribed annulus in cycles per meter

% Define coordinates with zeropad
k_v = center_freq + ((1./ss(2)) * [-1 1] / 2);
k_v = linspace(k_v(1), k_v(2), size(complex_data,2));
k_u = (1./ss(1)) * [-1 1] / 2;
k_u = linspace(k_u(1), k_u(2), size(complex_data,1));

% Coordinates onto which we will interpolate
new_data_size = size(complex_data); % Assume output array will have same size as input array
                                    % This datasize assumption is not
                                    % necessary, especially for significant
                                    % inscription loss, but it should
                                    % result in sufficient sample density
                                    % for all cases.
new_r = linspace(k_r(1), k_r(2), new_data_size(2)).';
new_a = linspace(k_a(1), k_a(2), new_data_size(1)).';

%% Resample
% We do the interpolation in two steps.  First we interpolate across
% columns of contant U to constant polar angle-- rectangular to
% "keystone"), and then we interpolate across rows of constant V to
% constant frequency-- "keystone" to polar).
keystone = zeros(numel(new_a), size(complex_data,2));
for row=1:size(complex_data,2)
    new_u = tan(new_a) * k_v(row);
    keystone(:,row) = interp1(k_u, complex_data(:,row), new_u, 'linear');
end
keystone = keystone.'; % We always want to operate along the first dimension
resampled_data = zeros(numel(new_r),numel(new_a));
for angle=1:size(keystone,2)
    k_r_keystone = k_v./cos(new_a(angle));
    resampled_data(:,angle) = interp1(k_r_keystone, keystone(:,angle), new_r, 'linear');
end
resampled_data = resampled_data.'; % Return to original orientation

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////