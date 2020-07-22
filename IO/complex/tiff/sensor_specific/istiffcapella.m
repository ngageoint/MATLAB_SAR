function [ boolout ] = istiffcapella( filename )
%ISTIFFCAPELLA Test if data is Capella TIFF format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


tags = read_tiff_tags(filename);
boolout = isfield(tags{1}, 'Software') && strncmpi(tags{1}.Software,'Capella',7);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////