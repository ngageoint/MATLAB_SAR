function [ output_data ] = deweightmem( input_data, weight_fun, oversample_rate, dim )
%DEWEIGHTMEM Make complex SAR uniformly weighted in one dimension
%
%    output_data = deweightmem(input_data, weight_fun, oversample_rate, dim)
%
%       Parameter name    Description
% 
%       input_data        Array of complex values to deweight
%       weight_fun        Description of weighting applied.  Either a
%                            function handle to a function that takes a
%                            single argument (number of elements) and
%                            produces the weighting to apply, or a vector
%                            that is the weighting function sampled.
%       oversample_rate   Amount of sampling beyond the ImpRespBW in the
%                            processing dimension. (Default is Nyquist
%                            sampling = 1).
%       dim               Dimension over which to perform deweighting.
%                            Default is 1.
%       output_data       INPUT_DATA with normalization applied
%
% This implementation assumes that the data has already been "deskewed" and
% that the frequency support is centered.
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('oversample_rate','var')
    oversample_rate = 1;
end
if ~exist('dim','var')
    dim = 1;
end

data_size = size(input_data,dim);
weight_size = round(data_size/oversample_rate); % Weighting only valid across ImpRespBW
if ~exist('weight_fun','var')||isempty(weight_fun) % No weighting passed in.  Do nothing.
    output_data = input_data;
    return;
elseif isa(weight_fun, 'function_handle')
    weighting = weight_fun(weight_size);
elseif isvector(weight_fun)
    weighting = interpft(weight_fun,weight_size);
end
% weighting = weighting * sqrt(mean(1./weighting.^2)); % We want total power maintained
weight_zp = ones(data_size,1); % Don't scale outside of ImpRespBW
weight_zp(floor((data_size-weight_size)/2)+(1:weight_size)) = weighting;
if dim==2
    weight_zp = weight_zp.';
end

% Divide out weighting in spatial frequency domain
output_data = fftshift(fft(input_data,[],dim),dim);
output_data = bsxfun(@rdivide,output_data,weight_zp);
output_data = ifft(ifftshift(output_data,dim),[],dim);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////