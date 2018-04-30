function target = camera(this,long,lat,alt,heading,tilt,roll,varargin)
%KML.camera(long,lat,alt) Write a Camera position KML with lat, long, alt,
% heading, tilt, roll
%  
%   Based on work by Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   Original Work: Copyright 2012 Rafael Fernandes de Oliveira
%   Author: Thomas Fell (t.r.fell@liv.ac.uk)
%   Version: 1.0 $  $Date: 2013/04/18 20:00:00 $
  
    target = struct('type','','id','');

    [long,lat,heading,tilt,roll] = this.checkUnit(long,lat,heading,tilt,roll);

    p = inputParser;
    
    nlat = numel(lat);
    
    p.addRequired('lat',@(a)isnumeric(a) && isvector(a) &&~isempty(a));
    p.addRequired('long',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('alt',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('heading',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('tilt',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('roll',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('id',kml.getTempID('kml_camera'),@ischar);
    p.addParamValue('name','kml_camera',@ischar);
    p.addParamValue('altitudeMode','absolute',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));
    
    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);
    
    p.parse(lat,long,alt,heading,tilt,roll,varargin{:});
    
    arg = p.Results;
    
    
    placemark   = this.xml.createElement('Placemark');
    camera       = this.xml.createElement('Camera');
    
    placemark.setAttribute('id',arg.id);
    placemark.appendChild(this.textNode('name',arg.name));
    
    if ~isempty(arg.timeStamp)
        timeStamp = this.xml.createElement('TimeStamp');
        timeStamp.appendChild(this.textNode('when',arg.timeStamp));
        placemark.appendChild(timeStamp);
    elseif ~isempty(arg.timeSpanBegin) || ~isempty(arg.timeSpanEnd)
        timeSpan = this.xml.createElement('TimeSpan');
        if ~isempty(arg.timeSpanBegin)
            timeSpan.appendChild(this.textNode('begin',arg.timeSpanBegin));
        end
        
        if ~isempty(arg.timeSpanEnd)
            timeSpan.appendChild(this.textNode('end',arg.timeSpanEnd));
        end
        placemark.appendChild(timeSpan);
    end
    
    camera.appendChild(this.textNode('latitude',num2str(arg.lat)));
    camera.appendChild(this.textNode('longitude',num2str(arg.long)));
    camera.appendChild(this.textNode('altitude',num2str(arg.alt)));
    camera.appendChild(this.textNode('heading',num2str(arg.heading)));
    camera.appendChild(this.textNode('tilt',num2str(arg.tilt)));
    camera.appendChild(this.textNode('roll',num2str(arg.roll)));
    camera.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
    
    placemark.appendChild(camera);
    this.doc.appendChild(placemark);
    
    target.type = 'Placemark';
    target.id   = arg.id;
end