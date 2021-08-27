function boolout = iscrsd( filename )
%ISCRSD Checks to see whether input file is in the CRSD format
%
% Written by: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b','UTF-8');
boolout = strncmp(fgets( fid, 20 )','CRSD/1.0.0',8);
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////