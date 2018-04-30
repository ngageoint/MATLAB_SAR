function subaperturedemo(varargin)
%SUBAPERTUREDEMO Read image, generate, and view subaperture-processed image using MATLAB GUIs.
%    subaperturedemo('PropertyName',PropertyValue,...)
%
% Calculates the subaperture-processed image of a complex image (selected
% through a MATLAB dialog box) using the properties specified. The AOI over
% which to compute the subaperture-processed image is selected through a MATLAB
% GUI. This version of subaperture processing does NOT require that the complete
% data fit into memory. It processes from any format handled by OPEN_READER and
% outputs to files in SIO format.
%
%       Property name     Description
%       frames            number of frames (default = 7)
%       apfraction        fraction of aperture for each subaperture
%                            (default = .25)
%       method            'normal' (default), 'fullpixel', or 'minimal'
%       platformdir       platform direction, 'right' (default) or 'left'
%       dim               dimension over which to split subaperture
%                            (default = 1)
%       fill              fill factor
%
% Output frames are stored in one frame per file, where the first frame has
% the filename OUTFILE, and each consecutive frame has the frame number
% appended onto the filename.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

demo_core(@subaperturefile, varargin{:});

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////