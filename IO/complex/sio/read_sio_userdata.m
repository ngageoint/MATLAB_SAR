function [ userdata_struct, userdata_length_in_bytes ] = read_sio_userdata( filename, endian )
%READ_SIO_USERDATA Extracts user data from SIO files
%
% Assumes that you already know the endianness ('b' for big-endian, 'l' for
% little) and that this file has user data
%
% Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r',endian,'UTF-8');
ihdr = fread(fid,6,'*uint32'); % Read SIO header
userdata_length_in_bytes=4; % uint32 telling how many pairs of user data
userdata_struct = struct();
for i=1:ihdr(6)
    namebytes = fread(fid,1,'*uint32'); % 4 bytes
    name=fread(fid,namebytes,'*char')';
    
    valuebytes = fread(fid,1,'*uint32'); % 4 bytes
    value=fread(fid,valuebytes,'*char')';
    
    userdata_struct.(genvarname(name))=value;
    userdata_length_in_bytes=userdata_length_in_bytes+namebytes+valuebytes+8;
end
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////