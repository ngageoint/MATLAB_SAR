function [ boolout ] = isgff( filename )
%ISGFF Sandia GSAT Image File Format
%
% Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b'); 
format=fread(fid,7,'uchar=>char')';
fclose(fid);
boolout = strcmp(format,'GSATIMG');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////