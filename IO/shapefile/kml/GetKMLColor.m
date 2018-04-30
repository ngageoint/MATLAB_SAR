function Color = GetKMLColor(InColor,Transparency)
%GETKMLCOLOR Convert color description into KML format
%
%Color in kml is defined as TTBBGGRR, where TT is transparency.  All values
%are 0-255 in hex
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if length(InColor) == 1
    Color = '00000000';
    return;
end

%convert transparency to 00-ff in Hex
TransInt = floor(Transparency*255);
TransStr = sprintf('%2.2x',TransInt);

%if RGB Triplet is passed in, convert to hex and build color string
if (InColor(1) >= 0 && InColor(1) <= 1 &&...
     InColor(2) >= 0 && InColor(2) <= 1 &&...
     InColor(3) >= 0 && InColor(3) <= 1)
    R = sprintf('%2.2x',floor(InColor(1)*255));
    G = sprintf('%2.2x',floor(InColor(2)*255));
    B = sprintf('%2.2x',floor(InColor(3)*255));
    Color= strcat(TransStr,B,G,R);
    return;
end

switch lower(InColor)
    case 'red'
        ColorString = '0000ff';
    case 'orange'
        ColorString = '0066ff';
    case 'yellow'
        ColorString = '00ffff';
    case 'green'
        ColorString = '00ff00';
    case 'cyan'
        ColorString = 'ffff00';
    case 'blue'
        ColorString = 'ff0000';
    case 'magenta'
        ColorString = 'ff00ff';
    case 'white'
        ColorString = 'ffffff';
    case 'black'
        ColorString = '000000';
end

Color = strcat(TransStr,ColorString);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////