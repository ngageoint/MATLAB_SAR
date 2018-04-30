function boolout = iscphd30( filename )
%ISCPHD30 Checks to see whether input file is in the CPHD 3.0 format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

boolout = false;
if exist(filename,'file')==2
    % Do quick easy check first
    fid = fopen(filename,'r','b');
    boolout = strncmp(fread( fid, 11, '*char' )','BegPreamble',11);
    fclose(fid);
    % Then check more thoroughly
    if boolout
        try
            preamble = read_cphd_preamble(filename);
            boolout = preamble.Version >= 2; % Version 2 seems to be similar enough to 3 to work for most things
        catch
            boolout = false;
        end
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////