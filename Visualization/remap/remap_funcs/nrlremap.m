function out = nrlremap(x, a, c)
% NRLREMAP Lin-log style remap
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin<3
    c = 220; % Output "knee" in lin-log curve
end
if nargin<2
    a = 1; % Scale factor of 99th percentile for input "knee"
end

A = abs(single(x));
Amin = min(A(:));
P99 = prctile(A(:),99);
b = (255-c)/log10((max(A(:))-Amin)/((a*P99)-Amin));

out = uint8(zeros(size(x)));
linear_region = (A<=a*P99);
out(linear_region) = (A(linear_region) - Amin)*c/((a*P99)-Amin);
out(~linear_region) = b*log10((A(~linear_region)-Amin)/((a*P99)-Amin))+c;

end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
