function [ output_struct ] = parsepdvinputs( varargin )
%PARSEPDVINPUTS Parse input arguments for pdv
%
%   Used in both in-memory and in-file version of PDV
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('deltax', 0.25, @(x) isscalar(x)&&(x>0));
p.addParamValue('filtersize', 5, @(x) isvector(x)&&(length(x)<=2)&&all(x>0));
p.addParamValue('filtertype', 'mean', @(x)any(strcmpi(x,{'mean','median'})));
p.addParamValue('dim', 1, @(x) isequal(x,1)||isequal(x,2));
p.addParamValue('framenumber',1);
p.FunctionName = 'PDV';
p.parse(varargin{:});

output_struct=p.Results;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////