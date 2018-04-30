function [AzZeroPad, RnZeroPad] = ComputeZeroPad(Data)
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%Az Zero Pad
MeanData = mean(abs(Data));
%look for ramp up and ramp down
AzZeroPad = GetZeroPad(MeanData);

%Range Zero Pad
MeanData = mean(abs(Data),2);
%look for ramp up and ramp down
RnZeroPad = GetZeroPad(MeanData);


function ZeroPad = GetZeroPad(data)

%some statistics
nx = length(data);
diffeddata = abs(diff(data));

meanval = mean(diffeddata(:));
stdval = std(diffeddata(:));

StartVal = 0;
%start at beginning and look for large jump
for i=1:nx
  if (diffeddata(i) > meanval)
    StartVal = i;
    break;
  end
  if (diffeddata(i) > meanval+stdval)
    StartVal = i;
    break;
  end
end

if (StartVal == 0)
  StartVal = 1;
end

ZeroPad = nx/(nx-2*StartVal);

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////