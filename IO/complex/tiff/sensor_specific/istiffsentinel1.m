function [ boolout ] = istiffsentinel1( filename )
%ISTIFFSENTINEL1 Test if data is Sentinel-1 TIFF
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

pathname = fileparts(filename);
s1xmlfile = dir(fullfile(pathname,'..','manifest.safe')); % Sentinel 1 manifest file
boolout = (length(s1xmlfile)==1)&&issentinel1slc(fullfile(fileparts(pathname),s1xmlfile.name));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////