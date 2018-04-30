function boolout = istiff(filename)
% Tagged image format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename,'r'); 
endianflag = fread(fid,2,'uchar=>char')'; % Big or little endian???
if strcmp(endianflag,'II')
    endian='l';
elseif strcmp(endianflag,'MM')
    endian='b';
else
    boolout=false;
    fclose(fid);
    return;
end
magicnumber = fread(fid,1,'int16',0,endian);
fclose(fid);
boolout = (magicnumber==42);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////