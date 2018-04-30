function [ pdvout ] = pdvmem( compleximage, varargin )
%PDVMEM Phase Derivative Value
%    pdvout = pdv(compleximage, 'PropertName',PropertyValue,...) calculates
%    the phase derivate of compleximage stored in memory using the
%    property specified:
%
%       Property name     Description
%       deltax            pixel shift (default = 0.25)
%       filtersize        size of smoothing filter (default = [5 5])
%       filtertype        type of filter, 'mean' (default) or 'median'
%       dim               dimension over which to calculate phase gradient
%                            (default = 1)
%
% Written by: Wade Schwartzkopf
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse and validate arguments if not already done
if isreal(compleximage) || (~isnumeric(compleximage)) || (ndims(compleximage) > 2)
    error('Input image must be a single complex image');
end
if (nargin>1)&&isstruct(varargin{1})
    inargs=varargin{1};
else
    inargs=parsepdvinputs(varargin{:});
end

%% Setup data arrays
prlength=size(compleximage,inargs.dim);
phaseramp=fftshift(exp(1i*pi*inargs.deltax*(-prlength/2:prlength/2-1)/prlength));
if inargs.dim==1 % Process vertically
    phaseramp=phaseramp.';
    phaserampmatrix=repmat(phaseramp,1,size(compleximage,2));
else % horizontally
    phaserampmatrix=repmat(phaseramp,size(compleximage,1),1);
end

%% Calculate phase gradient
cifft=fft(compleximage,[],inargs.dim);
shift1=ifft(cifft.*conj(phaserampmatrix),[],inargs.dim);
shift2=ifft(cifft.*phaserampmatrix,[],inargs.dim);
clear cifft phaserampmatrix;
pdvout=shift1.*conj(shift2); % phase derivative
clear shift1 shift2

%% "Boxcar" filtering
if any(inargs.filtersize>1)
    if(strcmpi(inargs.filtertype,'mean')) % Mean filter
        filtercoefs=ones(inargs.filtersize); filtercoefs=filtercoefs/numel(filtercoefs);
%         pdvout=conv2(padarray(pdvout,fix(size(filtercoefs)/2),'replicate','both'),filtercoefs,'valid');
        pdvout=conv2(pdvout,filtercoefs,'same'); % Assumes zeros outside edge of image
    else % Median filter (NRL appears to median filter angle, not real and complex parts)
        pdvout=complex(medfilt2(real(pdvout),inargs.filtersize.*[1 1]),...
            medfilt2(imag(pdvout),inargs.filtersize.*[1 1]));
    end
end

pdvout=angle(pdvout)/inargs.deltax;
% pdvout((abs(real(pdvfiltered))<.0000001)&(abs(imag(pdvfiltered))<.0000001))=0; % Remove spurious phases

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////