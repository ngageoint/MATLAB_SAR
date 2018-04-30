function fwriten(fid,a,N)
%FWRITEN Write a value to a file in text in a fixed length field.
% 
% If A is an integer, it is padded with zeros to the left.  If A is a char,
% it is padded with spaces to the right.  N is the length (in characters)
% of the field.
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if isa(a,'integer')
    % This is the reasonable way to zeropad integers using sprintf, but it
    % doesn't seem to work for very large values of 'a' for versions of
    % MATLAB prior to 2011b.  Must have been a bug in the Matlab sprintf
    % implementation.
    % x = sprintf('%0*d', N, a);
    % The following workaround seems to work for all sizes of integers even
    % on versions of MATLAB prior to the 2011b fix.
    x = sprintf('%0*.0f', N, a);
else
    x = sprintf('%-*s', N, a);
end
fwrite(fid,x(1:N),'uint8');
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////