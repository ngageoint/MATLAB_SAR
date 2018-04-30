function f = newAnimation(this,animName)
%KML.NEWANIMATION(folderName) Creates an animation storyboard inside an kml (or another folder).
%  Example of use:
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    if nargin < 2
        animName = 'Unnamed Animation';
    end
    f = kmlAnimation(this,animName);
end