function [ boolout ] = istiffrs2( filename )
%ISTIFFRS2 Test if data is RADARSAT-2 TIFF
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% RS2 TIFF doesn't have much metadata in the TIFF itself, but it should
% have a co-located product.xml file
pathname = fileparts(filename);
rs2xmlfile = dir(fullfile(pathname,'product.xml')); % RADARSAT-2 metadata file
boolout = (length(rs2xmlfile)==1)&&isrs(fullfile(pathname,rs2xmlfile.name));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////