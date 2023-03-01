function [ boolout ] = iscsm( filename )
%ISCSM Is Cosmo Skymed HFD5 file format
%
% Also handles KOMPSAT-5, which is nearly identical format
%
% Written by: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b');
boolout = isequal(fread( fid, 8, 'uint8' )',[137 72 68 70 13 10 26 10]); % HDF5 signature
fclose(fid);
% Now that we know it is HDF5, is it CSM HDF5?
if boolout
    try % Error if root attribute 'Mission ID' does not exist
        boolout=any(strncmpi(h5readatt(filename,'/','Mission ID'), {'CSK','CSG','KMP'}, 3));
    catch % h5readatt replaced hdf5read at some point in MATLAB history
        try
            boolout=any(strncmpi(hdf5read(filename,'/','Mission ID'), {'CSK','CSG','KMP'}, 3));
        catch % Is there a more direct way to check for attribute existence first?
            boolout=false;
        end
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////