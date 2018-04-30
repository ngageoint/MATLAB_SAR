function [ boolout ] = iscos( filename )
%ISCOS TerraSAR-X COSAR file format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b');
fseek(fid,28,'bof');
boolout = strncmp( fread( fid, 4, 'uchar=>char' ), 'CSAR', 4);
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////