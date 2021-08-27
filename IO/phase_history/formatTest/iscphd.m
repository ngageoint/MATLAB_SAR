function boolout = iscphd( filename )
%ISCPHD Checks to see whether input file is in the CPHD format
%
% Written by: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b','UTF-8');
[cphd_str, version_str] = strtok(fgets( fid, 20 ),'/'); % First line is just CPHD/version
vers_parts = str2double(regexp(version_str(2:end),'\.','split'));
boolout = strcmp(cphd_str,'CPHD') && ...  % CPHD"X"
    ((vers_parts(1)==0&&vers_parts(2)>=3)||vers_parts(1)>0);  % 0.3 and above
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////