function [ im0_RGB ] = csimem(im0,varargin) 
% CSIMEM Color Subaperture Image
%    im0_RGB = csimem(compleximage, 'PropertyName', PropertyValue, ...)
%
% Displays subaperture information as color on full resolution data.
%
%       Property name     Description
%       dim               dimension over which to split subaperture
%                            (default = 1)
%       fill              fill factor (default = 1)
%       platformdir       platform direction, 'right' (default) or 'left'.
%                            Assumption is that 2nd dimension is increasing
%                            range.
%       fullres           boolean.  Determines whether to reapply original
%                            intensity to final color output.  This
%                            maintains full resolution from original input
%                            image.  (default = true)
%
% Authors: Chris Coleman, Mike Brennan
% Adapted for the MATLAB SAR Toolbox by Wade Schwartzkopf
%
% Assumptions:
% Input array im0 is a complex valued SAR map in the image domain.
%
% The output image is a 3-band complex image representing red, blue, and
% green.  For a basic remap that will convert the output to a form that the
% MATLAB viewers like imagesc can display, use densityremap, included in
% the MATLAB SAR toolbox.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('platformdir', 'right', @(x)any(strcmpi(x,{'right','left'})));
p.addParamValue('dim', 1, @(x) isequal(x,1)||isequal(x,2));
p.addParamValue('fill', 1, @(x) isscalar(x));
p.addParamValue('fullres', true, @(x) isscalar(x));
p.FunctionName = 'CSI';
p.parse(varargin{:});

if p.Results.fill>=1
    fill = p.Results.fill;
else
    fill = 1;
end

% CSI is always done along the first dimension 
% For range CSI, we simply temporarily swap array dimensions.
% This isn't very efficient but keeps the code simple
if p.Results.dim==2
    im0=im0.';
end

% Setup the 3 band of color filters
cmap = jet_wrapped(round(size(im0,1)/fill)); % jet-like colormap, used for filtering

% Move to the phase history domain
ph0 = fftshift(ifft(ifftshift(im0,1)),1); %Take the inverse FFT with shifting

% Apply the subap filters
ph_indices=floor((size(im0,1)-size(cmap,1))/2)+(1:size(cmap,1));
ph0_RGB=zeros(length(ph_indices), size(im0,2), 3);
ph0_RGB(:,:,1) = ph0(ph_indices,:) .* repmat(cmap(:, 1),1,size(im0,2)); % Red
ph0_RGB(:,:,2) = ph0(ph_indices,:) .* repmat(cmap(:, 2),1,size(im0,2)); % Green
ph0_RGB(:,:,3) = ph0(ph_indices,:) .* repmat(cmap(:, 3),1,size(im0,2)); % Blue

% Shift phase history to avoid having zeropad in middle of filter.  This
% fixes the purple sidelobe artifact.
filtershift=ceil(size(im0,1)/(4*fill));
ph0_RGB(:,:,1) = circshift(ph0_RGB(:,:,1), -filtershift); % Red
% Green already centered
ph0_RGB(:,:,3) = circshift(ph0_RGB(:,:,3), filtershift); % Blue

% FFT back to the image domain
im0_RGB = ifftshift(fft(fftshift(ph0_RGB,1),size(im0,1)),1);

% Reorient images as necessary
if p.Results.dim==2
    im0_RGB=permute(im0_RGB,[2 1 3]);
    im0=im0.';
elseif strcmp(p.Results.platformdir,'right')
    temp=im0_RGB(:,:,1);
    im0_RGB(:,:,1)=im0_RGB(:,:,3);
    im0_RGB(:,:,3)=temp;
end

% We could stop here, but some worry that we have lost resolution with the
% color filters.  We will replace the intensity with the original image
% intensity to main full resolution (in intensity, but not in color).
if p.Results.fullres
    % What we are doing conceptually
    % im0_HSV=rgb2hsv(abs(im0_RGB));
    % im0_HSV(:,:,3)=abs(im0); % Replace v with original intensity
    % % Remap could be done here as well
    % im0_RGB=hsv2rgb(im0_HSV);

    % An equivalent, but much more efficient way to do it
    scale_factor=abs(im0)./max(abs(im0_RGB),[],3);
    im0_RGB=bsxfun(@times,scale_factor,abs(im0_RGB));
end

end

% Copies functionality of jet colormap.  However, "wraps" around blue and
% red filters, so that filters are identically shaped (in a periodic
% sense).
function J = jet_wrapped(m)
n = ceil(m/4);
u = [(1:1:n)/n ones(1,n-1) (n:-1:1)/n]';
g = round((m-length(u))/2) + (1:length(u))';
r = g + n;
b = g - n;
g=mod(g-1,m)+1;
r=mod(r-1,m)+1;
b=mod(b-1,m)+1;
J = zeros(m,3);
J(r,1) = u(1:length(r));
J(g,2) = u(1:length(g));
J(b,3) = u(end-length(b)+1:end);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////