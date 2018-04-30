function RGBData = Pauli(HH,HV,VH,VV)

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

RedData = abs(HH - VV);
GreenData = abs(HV + VH);
BlueData = abs(HH + VV);

RedData = densityremap(RedData);
GreenData = densityremap(GreenData);
BlueData = densityremap(BlueData);

RGBData(:,:,1) = RedData;
RGBData(:,:,2) = GreenData;
RGBData(:,:,3) = BlueData;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////