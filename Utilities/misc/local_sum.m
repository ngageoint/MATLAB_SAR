function value = local_sum(input, window_size)
% LOCAL_SUM Fast computation of 2D local sum using integral images.
% 
% Algorithm is described in Viola and Jones, "Robust real-time face
% detection", International Journal of Computer Vision, 2004.  Computation
% time is independent of window size. Precision may be an issue for very
% large images, but a few thousand by a few thousand images seems to work
% with no issues.
%
% Suprisingly, for small window sizes, MATLAB's conv2 actually seems to be
% a faster way to do this, but this implementation gets much faster as
% window sizes get larger.  MATLAB Central's fastrunmean demonstrates
% another way of computing fast means/sums.  Fastrunmean is not quite as
% fast as this integral image implementation, but doesn't have the same
% potential precision issues for very large images.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

padsize = floor(window_size/2);
padded = padarray(input, padsize+1, 0, 'pre'); % Zeropad to handle local sums at edges
padded = padarray(padded, padsize, 0, 'post'); % Zeropad to handle local sums at edges
padded = cumsum(cumsum(double(padded),1),2); % Create integral image
value = padded(1:size(input,1),1:size(input,2))-... % Compute sum over local window
    padded(1:size(input,1),(window_size(2)+1):end)-...
    padded((window_size(1)+1):end,1:size(input,2))+...
    padded((window_size(1)+1):end,(window_size(2)+1):end);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////