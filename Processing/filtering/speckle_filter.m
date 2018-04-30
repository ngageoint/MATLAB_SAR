function [fil, b] = speckle_filter(I, varargin)
% Lee, Kuan, Frost, or MAP speckle filtering 
%    fil = speckle_filter(I, 'radius', radius, 'ENL', 1, 'filter', 'Lee')
%
% Calculates the Lee or Kuan speckle filter on the input image I.  I may be real or complex.  If real, it is
% assumed to be the amplitude (absolute value) of the image, without any remap applied.  If size(I,3) > 1, 
% then multiple coregistered input images are assumed, and each is filtered with the same filter coefficients.
%
%       Property name     Description
%       radius            size of the window to use for the filter.  The window is (2*radius+1) x (2*radius+1)
%                             (default is 2) 
%       ENL               Number of equivalent looks to use.  If not specified, it is estimated by assuming at
%                             least 5% of the image is homogenous.  Note: this is not used for the Frost filter.
%       filter            'Boxcar', 'Lee', 'Kuan' (default), 'MAP', or 'Frost'
%       frost_param       Positive integer giving the size of the Gaussian used in the Frost filter.  Ignored
%                             for all other filters.  Default is 1.
%
% Output is returned as a real matrix whose values are amplitude.  Bias is corrected before output.  
%
% Limitation: Will not perform as expected on remapped data
%
% Written by: Tom Braun, NGA/IBR, 703-735-2644 or Thomas.R.Braun@nga.mil, Thomas.R.Braun@nga.ic.gov
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Modified by Pat Cutler NGA (20140306): MAP computation

p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('radius', 2, @(x) isscalar(x)&&(x>=1));
p.addParamValue('ENL', 0, @(x) isscalar(x)&&(x>=0));
p.addParamValue('filter', 'Kuan', @(x) any(strcmpi(x,{'Boxcar', 'Lee', 'Kuan', 'Frost', 'MAP'})));
p.addParamValue('frost_param', 1, @(x) isscalar(x)&&(x>0));
p.FunctionName = 'speckle_filter';
p.parse(varargin{:});
radius = p.Results.radius;
ENL = p.Results.ENL;
frost_param = p.Results.frost_param;
lee = strcmpi(p.Results.filter, 'Lee');
kuan = strcmpi(p.Results.filter, 'Kuan');
frost = strcmpi(p.Results.filter, 'Frost');
map = strcmpi(p.Results.filter, 'MAP');
boxcar = strcmpi(p.Results.filter, 'Boxcar');

[x y z]=size(I);
I = double(abs(I)).^2;

if ~boxcar
    m = zeros(x,y,z);
    v = zeros(x,y,z);
    for i = -radius:radius
        for j = -radius:radius
            m(max(1,1+i):min(x, x+i), max(1,1+j):min(y,y+j),:) = m(max(1,1+i):min(x, x+i), max(1,1+j):min(y,y+j), :) + I(max(1,1-i):min(x,x-i), max(1,1-j):min(y,y-j),:);
        end
    end
    m = m ./ repmat(min(min([1:x]', ones(x,1)*(2*radius+1)), [x:-1:1]') * min(min([1:y], ones(1,y)*(2*radius+1)), [y:-1:1]), [1 1 z]);
    for i = -radius:radius
        for j = -radius:radius
            v(max(1,1+i):min(x, x+i), max(1,1+j):min(y,y+j), :) = v(max(1,1+i):min(x, x+i), max(1,1+j):min(y,y+j), :) + (I(max(1,1-i):min(x,x-i), max(1,1-j):min(y,y-j), :) - m(max(1,1-i):min(x,x-i), max(1,1-j):min(y,y-j), :)).^2;
        end
    end
    v = v ./ repmat(min(min([1:x]', ones(x,1)*(2*radius+1)), [x:-1:1]') * min(min([1:y], ones(1,y)*(2*radius+1)), [y:-1:1]), [1 1 z]);
    Ci2 = mean(v ./ (m.^2), 3);
end

if frost
    % small means potentially lead to large Ci2 values and cause a problem for the frost filter.  Eliminate them
    radius = ceil(-log(0.01)/(sqrt(2)*frost_param*min(Ci2(:)))); % make sure we average values until the exponential is below 0.01 at the given point
    fil = zeros(x,y,z);
    for k = 1:z
        fil(:,:,k) = c_frost(I(:,:,k), Ci2, radius, frost_param);
    end    
elseif ~boxcar   
    if ENL <= 0
        % This is an estimate that assumes that about 10% of the scene is homogeneous.  A much better estimate should probably be placed here
        tmp = Ci2(radius+1:end-radius,radius+1:end-radius); % skip the edges, since they have slightly different statistics
        tmp = sort(tmp(:));
        ENL = 1./mean(tmp(1:ceil(length(tmp).*0.1)));
    end
    Cn2 = 1./ENL;
    if lee
        Ct2 = Ci2 - Cn2;
    else
        Ct2 = (Ci2 - Cn2) ./ (1 + Cn2);
    end
    Ct2 = max(1e-20, Ct2);

    if lee || kuan
        fil = m + (I - m).*repmat(min(1, Ct2./Ci2), [1 1 z]);
    else % map
        %PJC 20140306 changes made according to Lopes et al. "Maximum A 
        %Posteriori Speckle Filtering and First Order Texture Models in SAR
        %Images" in Proc. IGARSS, College Park, MD. May 24-28, 1990. pp. 
        %2409-2412.
        %Gamma filtering is inly applicable to pixels within the conditions
        %specified in idx
        fil = I;
        idx = repmat(Ci2>=Cn2 & Ct2<=1/(1+ENL), [1 1 z]);
        Ct2 = repmat(Ct2,[1 1 z]);
        fil(idx) = (m(idx) .* 1./Ct2(idx) - ENL - 1 ...
            + sqrt((m(idx) .* 1./Ct2(idx) - ENL - 1).^2 ...
            + 4 .* ENL .* I(idx) .* m(idx) ./ Ct2(idx))) ...
            .* Ct2(idx) ./ 2;
    end
else
    fil = zeros(x,y,z);
    parfor k = 1:z
        fil(:,:,k) = conv2(I(:,:,k), ones(2*radius+1,2*radius+1)./((2*radius+1).^2), 'same');
    end
end

% convert back to amplitude and correct any bias issues
parfor k = 1:z
    fil(:,:,k) = fil(:,:,k).*(mean(mean(I(:,:,k)))./mean(mean(fil(:,:,k))));
end
fil = sqrt(fil);

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////