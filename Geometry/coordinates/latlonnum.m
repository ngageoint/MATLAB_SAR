function [ decimal_degrees ] = latlonnum( latlon_input )
%LATLONNUM Convert a variety of lat/long formats into decimal degree value
%
% Should handle any string compliant with the ISO 6709:2008 standard or any
% of a number of variants for describing lat/long coordinates.  Also
% handles degree/minutes/seconds passed in as a vector.
%
% Authors: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Check for empty string
if isempty(latlon_input)
    decimal_degrees = NaN;
    return;
end

%% Vector format degrees/minutes/seconds
if isnumeric(latlon_input)&&isvector(latlon_input)&&numel(latlon_input)<4
    decimal_degrees = sign(latlon_input(1))*polyval(abs(latlon_input(end:-1:1)),1/60);
    return;
end

%% String input
% Handles decimal degrees and degree/minutes/second with delimiters
delimiters = latlon_input(regexp(latlon_input,'[^.\d]')); % Any non-numeric characters in string
tokens_str = textscan(latlon_input,'%s','Delimiter',delimiters,'MultipleDelimsAsOne',true);
tokens = str2double(tokens_str{1});
decimal_degrees = polyval(abs(tokens(end:-1:1)),1/60);
if xor(any(latlon_input=='W'|latlon_input=='S'),any(latlon_input=='-'))
    decimal_degrees=-decimal_degrees;
end
% Handles degree/minutes/second with no delimiters DDD,DDDMM,DDDMMSS
if numel(tokens_str{1})==1
    for i=1:min(3,floor(numel(strtok(tokens_str{1}{1},'.'))/2)-1)
        decimal_degrees = fix(decimal_degrees/100)+rem(decimal_degrees,100)/60;
    end
end

% Error checking should occur here
if isempty(tokens)||numel(tokens)>3||...
    decimal_degrees<-180||decimal_degrees>360||...
    sum(isstrprop(latlon_input, 'alpha'))>1
    decimal_degrees = NaN; % Unparseable inputs are returned as NaN
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////