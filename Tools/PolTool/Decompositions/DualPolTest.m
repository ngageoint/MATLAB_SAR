function RGBData = DualPolTest(HH,HV,VH,VV)

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[ny nx] = size(HH);

%if only one transmit then use that....if quad we'll just use H-Tx
if (max(abs(HH(:))) > 0)
    RedTxAngle = 0;
    GreenTxAngle = 0;
    BlueTxAngle = 0;
    BlueRxAngle = 0;
else
    RedTxAngle = 90;
    GreenTxAngle = 90;
    BlueTxAngle = 90;
    BlueRxAngle = 90;
end

RedRxAngle = 135;
GreenRxAngle = 45;

%Red Channel
[HHCof HVCof VHCof VVCof] = ComputePolCoeff(RedTxAngle,RedRxAngle,0,0);
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
[HHCof HVCof VHCof VVCof] = ComputePolCoeff(GreenTxAngle,GreenRxAngle,0,0);
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
[HHCof HVCof VHCof VVCof] = ComputePolCoeff(BlueTxAngle,BlueRxAngle,0,0);
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

%scale data to mean + 3*sigma 
RedCutoff = mean(RedData(:)) + 3*std(RedData(:));
GreenCutoff = mean(GreenData(:)) + 3*std(GreenData(:));
BlueCutoff = mean(BlueData(:)) + 3*std(BlueData(:));

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