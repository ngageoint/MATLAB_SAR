function boolout = issicd01(filename)
% Sensor independent complex data (version 0.1)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b'); 
format=fread(fid,10,'uchar=>char')';
fclose(fid);
boolout = strcmp(format,'SICD-E/0.1');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////