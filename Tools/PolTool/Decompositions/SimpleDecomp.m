function RGBData = SimpleDecomp(HH,HV,VH,VV)

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

RedData = abs((HV+VH)/2);
GreenData = abs(HH);
BlueData = abs(VV);

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