function x4win = x4win(len)
%X4WIN Returns a 1/x^4 window of length x

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


x = linspace(-100,100,len);
y = x.^4;
y = y./max(y);
x4win = (1-y)';

end

