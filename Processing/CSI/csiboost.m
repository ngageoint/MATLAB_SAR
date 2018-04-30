function [ y ] = csiboost( a, boostlevel )
%CSIBOOST Optional remap for CSI that boost strongly colored pixels
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if (nargin<2)
    boostlevel=1; % 0 = no boost
end

%Detect image
y = abs(a);

hsvtest = rgb2hsv(y);
hsvmag = hsvtest(:,:,3);

%multiply Value by Saturation
hsvmag = hsvmag .* (hsvtest(:,:,2).^boostlevel);

%perform remap
stdmat = std(hsvmag(:));
meanmat = mean(hsvmag(:));
hsvmag(hsvmag > (meanmat+3*stdmat)) = meanmat+3*stdmat;
hsvmag = hsvmag./max(hsvmag(:));
hsvtest(:,:,3) = hsvmag;

%convert back to RGB
y = hsv2rgb(hsvtest);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////