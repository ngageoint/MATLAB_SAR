function [ ihdr, endian, data_offset, user_data ] = read_sio_meta( filename )
%READ_SIO_META Parse SIO header
%
% Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r','b');
ihdr = fread(fid,5,'*uint32');
fclose(fid);

switch ihdr(1)
    case hex2dec('FF017FFE') % Big endian, no user-data
        data_offset=20;
        endian='b';
        user_data=struct();
    case hex2dec('FE7F01FF') % Little endian, no user-data
        data_offset=20;
        endian='l';
        ihdr=swapbytes(ihdr);
        user_data=struct();
    case hex2dec('FF027FFD') % Big endian, with user-data
        endian='b';
        [user_data, user_data_length]=read_sio_userdata(filename,'b');
        data_offset=20+user_data_length;
    case hex2dec('FD7F02FF') % Little endian, with user-data
        endian='l';
        ihdr=swapbytes(ihdr);
        [user_data, user_data_length]=read_sio_userdata(filename,'l');
        data_offset=20+user_data_length;
    otherwise % Not an SIO file
        ihdr=[];
        endian=[];
        data_offset=0;
        user_data=struct([]);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////