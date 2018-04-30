function latlon_string = latlonstr(value, latlon, varargin)
% LATLONSTR Convert latitude/longitude decimal degree value to customizable string format
%
% LATLONSTR(DECIMAL_DEGREES, LATLON, 'PropertyName', PropertyValue)
% LATLONSTR([DEGREES MINUTES SECONDS], LATLON, 'PropertyName', PropertyValue)
%
% INPUTS:
%    VALUE:              Value of latitude or longitude in degrees or dms vector.
%    LATLON:             Either 'lat' or 'lon'.
%    Property name       Description
%       num_units        1 (decimal degrees), 2 degrees/minutes, 3
%                           degrees/minutes/seconds.   Default is 3.
%       delimiter        Single string or cell array of strings of
%                           separators between degrees/minutes/seconds/
%                           hemisphere.  Default is '' (empty).
%       include_symbols  True/false.  Whether to include °,'," symbols.
%                           Default is true.
%       signed           True/false.  Whether to use +/- or N/S/E/W to
%                           represent hemisphere.  Default is false
%                           (N/S/E/W).
%       precision        Number of decimal points shown in finest unit.
%                           Default is 5 if num_units==1, otherwise 0.
%       padded           True/false.  Whether to use zeros to pad out to
%                           consistent string length (3 digits for
%                           longitude degrees, 2 digits for all other
%                           elements).  Default is true.
%
% Supports ISO 6709:2008 formatted geographic coordinates
%    Annex D (human interface)
%       delimiter = ''; include_symbols = true; padded = true; signed = false
%    Annex H (string representation)
%       delimiter = ''; include_symbols = false; padded = true; signed = true
%
% Authors: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input args
if numel(value)>1 % Value passed in as degree/minutes/seconds
    value = latlonnum(value); % Assure we start with decimal degrees
end
p = inputParser; % Extract parameter-value pairs
p.addParamValue('num_units',3);
p.addParamValue('delimiter','');
p.addParamValue('include_symbols',true);
p.addParamValue('signed',false);
p.addParamValue('precision',[]);
p.addParamValue('padded',true);
p.FunctionName = mfilename;
p.parse(varargin{:});
% Precision.  Default is dependent on other input arguments.
if isempty(p.Results.precision)
    if p.Results.num_units==1
        precision = 5;
    else
        precision = 0;
    end
else
    precision = p.Results.precision;
end
% Symbols
if p.Results.include_symbols
    latlon_symbols = {char(176) '''' '"'};
else
    latlon_symbols = {'' '' ''};
end
% Delimiters
if iscell(p.Results.delimiter)
    delimiter = p.Results.delimiter;
else
    delimiter = repmat({p.Results.delimiter},[p.Results.num_units 1]);
end
if p.Results.signed
    delimiter{p.Results.num_units} = ''; % No separator needed for hemisphere
end

%% Differences between latitude and longitude
if strncmpi(latlon, 'lat', 3)
    if value(1)>0
        hemisphere = 'N';
        latlon_sign = '+';
    else
        hemisphere = 'S';
        latlon_sign = '-';
    end
    degrees_digits = 2;
elseif strncmpi(latlon, 'lon', 3)
    if value(1)>180
        value(1) = value(1) - 360;
    end
    if value(1)>0
        hemisphere = 'E';
        latlon_sign = '+';
    else
        hemisphere = 'W';
        latlon_sign = '-';
    end
    degrees_digits = 3;
end

%% Computer degree/minutes/seconds
new_value = abs(value);
for i = 1:p.Results.num_units
    fraction = rem(new_value,1);
    value(i) = fix(new_value);
    new_value = fraction*60;
end
value(end) = value(end) + fraction;
if p.Results.num_units > 1 && round(value(end),precision) == 60  % 60 seconds is invalid
    value(end) = 0;
    value(end-1) = value(end-1) + 1;
    if p.Results.num_units == 3 && value(end-1) == 60  % if adding 1 makes minutes 60
        value(end-1) = 0;
        value(1) = value(1) + 1;
    end
end

%% Build string
latlon_string = '';
for i=1:p.Results.num_units
    if p.Results.padded
        if i==1
            int_digits = degrees_digits;
        else
            int_digits = 2; % True for all but longitude degrees
        end
    else
        int_digits = 1;
    end
    if i==p.Results.num_units
        precision_digits = precision;
    else
        precision_digits = 0;
    end
    if precision_digits>0
        int_digits = int_digits + 1; % Account for the decimal point
    end
    latlon_string = sprintf('%s%0*.*f%s%s%s',latlon_string,...
        int_digits+precision_digits,precision_digits,abs(value(i)),...
        latlon_symbols{i},delimiter{i});
end
% Add hemisphere
if p.Results.signed
    latlon_string = [latlon_sign latlon_string];
else
    latlon_string = [latlon_string hemisphere];
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////