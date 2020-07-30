function [ boolout ] = isiceye( filename )
%ISCSM Is ICEYE HFD5 file format
%
% Written by: Tim Cox, NRL
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b');
boolout = isequal(fread( fid, 8, 'uint8' )',[137 72 68 70 13 10 26 10]); % HDF5 signature
fclose(fid);
% Now that we know it is HDF5, is it ICEYE?
if boolout
    try % Error if root attribute 'mission_name' does not exist
        boolout=strncmpi(h5read(filename,'/satellite_name'), 'ICEYE', 5) && ...
            strcmpi(h5read(filename,'/product_level'), 'SLC');
    catch % h5readatt replaced hdf5read at some point in MATLAB history
        try
            boolout=strncmpi(hdf5read(filename,'/satellite_name'), 'ICEYE', 5) && ...
                strcmpi(hdf5read(filename,'/product_level'), 'SLC');
        catch % Is there a more direct way to check for attribute existence first?
            boolout=false;
        end
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////