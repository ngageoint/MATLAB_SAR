function [ boolout ] = istiffrcm( filename )
%ISTIFFRCM Test if data is Radar Constellation Mission TIFF
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

pathname = fileparts(filename);
manifest_filename = dir(fullfile(pathname,'..','manifest.safe'));
product_filename = dir(fullfile(pathname,'..','metadata','product.xml'));
boolout = (length(manifest_filename)==1) && ...
    isrcmmanifest(fullfile(fileparts(pathname),manifest_filename.name)) && ...
    isrs(fullfile(fileparts(pathname),'metadata',product_filename.name));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////