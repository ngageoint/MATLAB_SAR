function RGBData = CreateDecomp(fname,HH,HV,VH,VV)

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%parse decomp data out of file
fid = fopen(fname,'r');

line = fgetl(fid);
index = strfind(line,'=');
RedTxAngle = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
RedTxEllip = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
RedRxAngle = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
RedRxEllip = str2double(line(index+1:end));
line = fgetl(fid);

line = fgetl(fid);
index = strfind(line,'=');
GreenTxAngle = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
GreenTxEllip = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
GreenRxAngle = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
GreenRxEllip = str2double(line(index+1:end));
line = fgetl(fid);

line = fgetl(fid);
index = strfind(line,'=');
BlueTxAngle = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
BlueTxEllip = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
BlueRxAngle = str2double(line(index+1:end));
line = fgetl(fid);
index = strfind(line,'=');
BlueRxEllip = str2double(line(index+1:end));

fclose(fid);

[ny nx] = size(HH);

%apply decomp to data

%Red Channel
[HHCof HVCof VHCof VVCof] = ComputePolCoeff(RedTxAngle,RedRxAngle,RedTxEllip,RedRxEllip);
im = complex(zeros(ny,nx),zeros(ny,nx));
if (abs(HHCof) > 0)
    im = im + HH*HHCof;
end
if (abs(HVCof) > 0)
    im = im + HV*HVCof;
end
if (abs(VHCof) > 0)
    im = im + VH*VHCof;
end
if (abs(VVCof) > 0)
    im = im + VV*VVCof;
end
RedData = abs(im);

%Green Channel
[HHCof HVCof VHCof VVCof] = ComputePolCoeff(GreenTxAngle,GreenRxAngle,GreenTxEllip,GreenRxEllip);
im = complex(zeros(ny,nx),zeros(ny,nx));
if (abs(HHCof) > 0)
    im = im + HH*HHCof;
end
if (abs(HVCof) > 0)
    im = im + HV*HVCof;
end
if (abs(VHCof) > 0)
    im = im + VH*VHCof;
end
if (abs(VVCof) > 0)
    im = im + VV*VVCof;
end
GreenData = abs(im);

%Blue Channel
[HHCof HVCof VHCof VVCof] = ComputePolCoeff(BlueTxAngle,BlueRxAngle,BlueTxEllip,BlueRxEllip);
im = complex(zeros(ny,nx),zeros(ny,nx));
if (abs(HHCof) > 0)
    im = im + HH*HHCof;
end
if (abs(HVCof) > 0)
    im = im + HV*HVCof;
end
if (abs(VHCof) > 0)
    im = im + VH*VHCof;
end
if (abs(VVCof) > 0)
    im = im + VV*VVCof;
end
BlueData = abs(im);

%scale data to mean + 2*sigma 
RedCutoff = mean(RedData(:)) + 2*std(RedData(:));
GreenCutoff = mean(GreenData(:)) + 2*std(GreenData(:));
BlueCutoff = mean(BlueData(:)) + 2*std(BlueData(:));

RedData(RedData > RedCutoff) = RedCutoff;
GreenData(GreenData > GreenCutoff) = GreenCutoff;
BlueData(BlueData > BlueCutoff) = BlueCutoff;

RedData = RedData./RedCutoff;
GreenData = GreenData./GreenCutoff;
BlueData = BlueData./BlueCutoff;

RGBData(:,:,1) = RedData;
RGBData(:,:,2) = GreenData;
RGBData(:,:,3) = BlueData;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////