function modelTour(this,time,long,lat,alt,heading,tilt,roll,varargin)
%THIS FUNCTION IS DEPRECATED - USE KML.NEWANIATION INSTEAD
%KML.MODELTOUR (long,lat,alt,heading,tilt,roll) Animates a 3D model (specified by the pair
%  attribute 'model','modelfile.dae') in the position and orientation given by long,
%  lat, alt, heading, tilt and roll. The time input controls how long the model takes
%  to fly from one coordinate to the next.
%
%  Many model files can be found in the <a href="http://sketchup.google.com/3dwarehouse/">Google 3D Warehouse</a> website.
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $
    
    [long,lat,heading,tilt,roll] = this.checkUnit(long,lat,heading,tilt,roll);
    
    p = inputParser;
    
    nlat = numel(lat);
    
    p.addRequired('time',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('lat',@(a)isnumeric(a) && isvector(a) &&~isempty(a));
    p.addRequired('long',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('alt',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('heading',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('tilt',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('roll',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('id','kml_modelTour',@ischar);
    p.addParamValue('name','kml_modelTour',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('model','',@ischar);
    p.addParamValue('scaleX',1,@(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('scaleY',1,@(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('scaleZ',1,@(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('scale',1, @(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));
    
    p.addParamValue('tourName','Play Me!',@ischar);
    p.addParamValue('cameraMode','above',@(a)ismember(a,{'behind','above','fixed'}));
    p.addParamValue('cameraDistance',1e3,@isnumeric);
    
    
    p.parse(time,lat,long,alt,heading,tilt,roll,varargin{:});
    
    arg = p.Results;
    
    
    if isempty(arg.model)
        error('Missing model parameter')
    end
    
    arg.scaleX = arg.scaleX .*arg.scale;
    arg.scaleY = arg.scaleY .*arg.scale;
    arg.scaleZ = arg.scaleZ .*arg.scale;
    
    if numel(arg.scaleX)==1 && nlat>1
        arg.scaleX = repmat(arg.scaleX,size(lat));
    end
    
    if numel(arg.scaleY)==1 && nlat>1
        arg.scaleY = repmat(arg.scaleY,size(lat));
    end
    
    if numel(arg.scaleZ)==1 && nlat>1
        arg.scaleZ = repmat(arg.scaleZ,size(lat));
    end
    
    
    locationID    = [arg.id '_location'];
    orientationID = [arg.id '_orientation'];
    scaleID       = [arg.id '_scale'];
    
    placemark   = this.xml.createElement('Placemark');
    model       = this.xml.createElement('Model');
    location    = this.xml.createElement('Location');
    orientation = this.xml.createElement('Orientation');
    scale       = this.xml.createElement('Scale');
    link        = this.xml.createElement('Link');
    
    placemark.setAttribute('id',arg.id);
    placemark.appendChild(this.textNode('name',arg.name));
    placemark.appendChild(this.textNode('visibility',num2str(arg.visibility)));
    placemark.appendChild(this.textNode('description',arg.description));
    
    location.setAttribute('id',locationID);
    location.appendChild(this.textNode('latitude',  num2str(lat(1),16)));
    location.appendChild(this.textNode('longitude', num2str(long(1),16)));
    location.appendChild(this.textNode('altitude',  num2str(alt(1),16)));
    
    orientation.setAttribute('id',orientationID);
    orientation.appendChild(this.textNode('heading', num2str(heading(1))));
    orientation.appendChild(this.textNode('tilt',    num2str(tilt(1))));
    orientation.appendChild(this.textNode('roll',    num2str(roll(1))));
    
    scale.setAttribute('id',scaleID);
    scale.appendChild(this.textNode('x', num2str(arg.scaleX(1))));
    scale.appendChild(this.textNode('y', num2str(arg.scaleY(1))));
    scale.appendChild(this.textNode('z', num2str(arg.scaleZ(1))));
    
    link.appendChild(this.textNode('href',arg.model));
    
    model.setAttribute('id',['Model_' arg.id]);
    model.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
    
    model.appendChild(location);
    model.appendChild(orientation);
    model.appendChild(scale);
    model.appendChild(link);
    
    placemark.appendChild(model);
    this.doc.appendChild(placemark);
    
    tour     = this.xml.createElement('gx:Tour');
    playlist = this.xml.createElement('gx:Playlist');
    tour.appendChild(this.textNode('name',arg.tourName));
    
    for i = 2:nlat-1
        
        dT = time(i)-time(i-1);
        
        % Update object position
        animUpdate  = this.xml.createElement('gx:AnimatedUpdate');
        update      = this.xml.createElement('Update');
        change      = this.xml.createElement('Change');
        location    = this.xml.createElement('Location');
        orientation = this.xml.createElement('Orientation');
        scale       = this.xml.createElement('Scale');
        
        animUpdate.appendChild(this.textNode('gx:duration',num2str(dT)));
        
        location.setAttribute('targetId',locationID);
        location.appendChild(this.textNode('latitude',  num2str(lat(i),16)));
        location.appendChild(this.textNode('longitude', num2str(long(i),16)));
        location.appendChild(this.textNode('altitude',  num2str(alt(i),16)));
        
        orientation.setAttribute('targetId',orientationID);
        orientation.appendChild(this.textNode('heading', num2str(heading(i))));
        orientation.appendChild(this.textNode('tilt',    num2str(tilt(i))));
        orientation.appendChild(this.textNode('roll',    num2str(roll(i))));
        
        scale.setAttribute('targetId',scaleID);
        scale.appendChild(this.textNode('x', num2str(arg.scaleX(i))));
        scale.appendChild(this.textNode('y', num2str(arg.scaleY(i))));
        scale.appendChild(this.textNode('z', num2str(arg.scaleZ(i))));
        
        change.appendChild(location);
        change.appendChild(orientation);
        change.appendChild(scale);
        
        update.appendChild(this.textNode('targetHref',''));
        update.appendChild(change);
        animUpdate.appendChild(update);
        
        playlist.appendChild(animUpdate);
        
        % Update camera position
        switch arg.cameraMode
            case 'fixed'
                wait  = this.xml.createElement('gx:Wait');
                wait.appendChild(this.textNode('gx:duration',num2str(dT)));
                playlist.appendChild(wait);
            case 'behind'
                dS = [long(i+1) - long(i-1); lat(i+1) - lat(i-1);];
                dS = arg.cameraDistance * dS./sqrt(sum(dS.^2))/1e6;
                
                flyTo = this.xml.createElement('gx:FlyTo');
                flyTo.appendChild(this.textNode('gx:duration',num2str(dT)));
                flyTo.appendChild(this.textNode('gx:flyToMode','smooth'));
                camera = this.xml.createElement('Camera');
                camera.appendChild(this.textNode('longitude',num2str(long(i)-dS(1),16)));
                camera.appendChild(this.textNode('latitude',num2str(lat(i)-dS(2),16)));
                camera.appendChild(this.textNode('altitude',num2str(alt(i),16)));
                h = heading(i)-180;
                
                if h < -360
                    h = 720+h;
                elseif h>360
                    h = h-720;
                end
                    
                camera.appendChild(this.textNode('heading',num2str(h)));
                camera.appendChild(this.textNode('tilt',num2str(90)));
                camera.appendChild(this.textNode('roll',num2str(0)));
                camera.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
                
                flyTo.appendChild(camera);
                playlist.appendChild(flyTo);
            case 'above'
                
                flyTo = this.xml.createElement('gx:FlyTo');
                flyTo.appendChild(this.textNode('gx:duration',num2str(dT)));
                flyTo.appendChild(this.textNode('gx:flyToMode','smooth'));
                lookAt = this.xml.createElement('LookAt');
                lookAt.appendChild(this.textNode('longitude',num2str(long(i),16)));
                lookAt.appendChild(this.textNode('latitude',num2str(lat(i),16)));
                lookAt.appendChild(this.textNode('altitude',num2str(alt(i) + arg.cameraDistance,16)));
                lookAt.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
                
                flyTo.appendChild(lookAt);
                playlist.appendChild(flyTo);
        end
    end
    
    tour.appendChild(playlist);
    this.doc.appendChild(tour);
    
    this.addIncludeFile(arg.model);
end