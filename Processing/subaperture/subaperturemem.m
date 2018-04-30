function result = subaperturemem(im, varargin)
%SUBAPERTUREMEM subaperture processed imaging
%    result = subaperturemem(compleximage, 'PropertyName', PropertyValue, ...)
%
% Calculates the subaperture-processed image on complex data that is held in memory
%
%       Property name     Description
%       frames            number of frames (default = 7)
%       apfraction        fraction of aperture for each subaperture
%                            (default = .25)
%       method            'normal' (default), 'fullpixel', or 'minimal'
%       platformdir       platform direction, 'right' (default) or 'left'
%       dim               dimension over which to split subaperture
%                            (default = 1)
%       fill              fill factor (default = 1)
%
% Output is stored in a cell array, where each element of the array is on
% frame.
%
% Limitation: Assumes complex data and all frames can be held in
% memory at once.
%
% Written by: Tom Braun
%             Wade Schwartzkopf
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse and validate arguments if not already done
if isreal(im) || (~isnumeric(im)) || (ndims(im) > 2)
    error('Input image must be a single complex image');
end
if (nargin>1)&&isstruct(varargin{1})
    inargs=varargin{1};
else
    inargs=parsesubapertureinputs(varargin{:});
end

%% Compute outpute size
num_x=size(im,inargs.dim); % Length along processing direction
if strcmp(inargs.method, 'minimal')
    inargs.output_res = ceil(num_x ./ inargs.fill ./ inargs.frames);
elseif strcmp(inargs.method, 'fullpixel')
    inargs.output_res = num_x;
else
    inargs.output_res = ceil(inargs.apfraction*num_x);
end

%% Setup for the subaperture processing
if (inargs.dim==2), im = im.'; end;
IM = fft(im);
clear im;
cutoff = floor(size(IM,1) ./ inargs.fill ./ 2);
IM = [IM(end-cutoff+1:end, :); IM(1:cutoff, :)]; % FFTSHIFT + chop off zeropad
result = cell(inargs.frames,1);
num_sa = ceil(inargs.apfraction * 2 * cutoff);
step = floor((2 * cutoff - num_sa)./(inargs.frames - 1));
offset = floor((2 * cutoff - (step*(inargs.frames-1)+num_sa)) ./ 2);
   
%% Run the subaperture processing.  We do it a frame at a time instead of a line at a time -- either works
for f = 1:inargs.frames
    if inargs.reverse_frames
        frame_num = inargs.frames - f + 1;
    else
        frame_num = f;
    end
    result{frame_num} = ifft(IM(offset + step*(f-1) + (1:num_sa), :), inargs.output_res);
    if inargs.dim == 2
        result{frame_num}=result{frame_num}.';
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////