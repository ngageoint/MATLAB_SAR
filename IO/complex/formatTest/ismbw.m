function boolout = ismbw(filename)
% Multi-band wrapper
%
% An intermediate format just for handling data that is contained in
% multiple files (color, polarimetric, etc.)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b');
id_string='MULTI-BAND WRAPPER';
firstline=fread(fid,length(id_string),'uchar=>char')';
boolout=strcmp(firstline,id_string);
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////