function boolout = iscnitf(filename)
% ISCNITF Check to see if a file is a complex NITF
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b'); 
FHDR = fread(fid,9,'uint8=>char')';
% fseek(fid,30,'cof');
% FTITLE = fread(fid,5,'uint8=>char')';
fclose(fid);
boolout=strcmp(FHDR,'NITF02.10')||strcmp(FHDR,'NITF02.00');
% No need to test for non-SICD NITF, since SICD test will be run first in
% guess_complex_format.m.
% boolout=boolout&&~strcmp(FTITLE,'SICD:');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////