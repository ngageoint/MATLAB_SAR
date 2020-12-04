function out = logremap(x)
% LOGREMAP Logarithmic remap
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

x = abs(single(x));
x = 10*x./mean(x(:)); % Fixes very small data (<< 1)
rcent  = 10*log10(sum(x(:).^2)/numel(x));
% Make minimum at least one
x = x + max(1 - min(x(:)),0);
x = 20*log10(x);

span_db = 50;
disp_min = max(min(x(:)),rcent - span_db/2);
disp_max = min(max(x(:)),rcent + span_db/2);

out = uint8(255*(x-disp_min)/(disp_max-disp_min));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////