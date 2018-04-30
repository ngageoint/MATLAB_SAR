function boolout = iscphdx( filename )
%ISCPHDX Checks to see whether input file is in the CPHD "X" format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b');
boolout = strncmp(fgets( fid, 20 )','CPHD/0.3',8);
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////