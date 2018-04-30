function [ acd_rgb ] = acdmem( reference_image, match_image, varargin )
%ACDMEM Amplitude change detection
%    acd_rgb = acdmem(reference_image, match_image, 'PropertyName', PropertyValue, ...)
%
% Compares two images by assigning each image to a different color and
% creating and RGB image of that composition.  Assumes both images are
% already registered with each other.  Input colors are in RGB space.
%
%       Property name     Description
%       reference_color   Color to which to assign reference image.
%                            Default = [1 0 0] (red).
%       match_color       Color to which to assign match image.  Default =
%                            [0 1 1] (cyan).
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('reference_color', [1 0 0], @isvalidrgb);
p.addParamValue('match_color', [0 1 1], @isvalidrgb);
p.FunctionName = 'acdmem';
p.parse(varargin{:});

rcolor=p.Results.reference_color;
mcolor=p.Results.match_color;

acd_rgb = cat(3, rcolor(1)*reference_image+mcolor(1)*match_image,...
                 rcolor(2)*reference_image+mcolor(2)*match_image,...
                 rcolor(3)*reference_image+mcolor(3)*match_image);

end

% RGB values are 3-element vectors between 0 and 1
function bool = isvalidrgb(in)
    bool=(max(size(in))==3)&&(min(size(in))==1)&&all(in>=0)&&all(in<=1);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////