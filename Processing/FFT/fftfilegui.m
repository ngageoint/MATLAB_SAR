function fftfilegui(varargin)
%FFTFILEGUI Read image, generate, and view data in the dispersed domain using MATLAB GUIs.
%
%    fftfile('PropertyName',PropertyValue,...)
%
% Function will query user for input file and AOI.  Only real parameter
% that can be passed in dimension for FFT to operate.  Default is 2D FFT.
%
%       Property name     Description
%       dims              dimensions along which to compute FFT (options:
%                            1, 2, or [1 2] (2-D FFT), default = [1 2])
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

demo_core(@fftfile, varargin{:});

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////