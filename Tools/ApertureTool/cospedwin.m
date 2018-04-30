function win = cospedwin(len,ped)
%COSPEDWIN Cosine on pedestal window 

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


%args: len - window length
%      ped - pedestal height (0-1)

if ped == 1 %this will blow up
    win = ones(len,1);
    return;
end

x = linspace(-pi,pi,len);
fac = 2./(1-ped);
win = ((cos(x)+1)./fac + ped)';

end

