function [ output_struct ] = parsesubapertureinputs( varargin )
%PARSESUBAPERTUREINPUTS Parse input arguments for subaperture processing
%
% Used in both in-memory and in-file version of subaperture processing
%
% Written by: Wade Schwartzkopf
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('frames', 7, @(x) isscalar(x)&&(x>=2)&&(mod(x,1)==0));
p.addParamValue('apfraction', 0.25, @(x) isscalar(x)&&(x<1)&&(x>0));
p.addParamValue('method', 'normal', @(x) any(strcmpi(x,{'minimal', 'normal','fullpixel'})));
p.addParamValue('platformdir', 'right', @(x)any(strcmpi(x,{'right','left'})));
p.addParamValue('dim', 1, @(x) isequal(x,1)||isequal(x,2));
p.addParamValue('fill', 1, @(x) isscalar(x)&&(x>=1));
p.FunctionName = 'SUBAPERTURE';
p.parse(varargin{:});

output_struct=p.Results;
output_struct.reverse_frames=(output_struct.dim==1)&&strcmp(output_struct.platformdir, 'right');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////