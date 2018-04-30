function PushpinFile = GetKMLPushpinFile(InColor)
%GETKMLPUSHPINFILE Get URL for pushpin based on color string
%
% Pushpins in Google Earth don't seem to use the 'color' tag, rather there
% are separate PNG files for each one.  This routine maps the color name
% passed in to the appropriate thumbtack file.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

switch lower(InColor)
    case 'red'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/red-pushpin.png';
    case 'yellow'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png';
    case 'green'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/grn-pushpin.png';
    case 'cyan'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/ltblu-pushpin.png';
    case 'blue'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/blue-pushpin.png';
    case 'magenta'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/pink-pushpin.png';
    case 'white'
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/wht-pushpin.png';
    case 'black'
        % No black thumbtack defined in Google Earth...
        PushpinFile = 'http://maps.google.com/mapfiles/kml/pushpin/wht-pushpin.png';
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////