function target = model(this,long,lat,alt,heading,tilt,roll,varargin)
%KML.MODEL(long,lat,alt,heading,tilt,roll) Places the 3D model (specified by the pair 
%  attribute 'model','modelfile.dae') in the position and orientation given by long,
%  lat, alt, heading, tilt and roll. If the inputs represent a vector of points, the 
%  model is repeated in each of them.
%  
%  Many model files can be found in the <a href="http://sketchup.google.com/3dwarehouse/">Google 3D Warehouse</a> website.
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

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
    
    p.addParamValue('id',kml.getTempID('kml_model'),@ischar);
    p.addParamValue('name','kml_model',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('model','',@ischar);
    p.addParamValue('scaleX',1,@(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('scaleY',1,@(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('scaleZ',1,@(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('scale',1, @(a)isnumeric(a) && (numel(a)==nlat || numel(a)==1));
    p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,alt,heading,tilt,roll,varargin{:});
    
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
    
    
    for i = 1:nlat
        target(i).type = 'Placemark';
        target(i).id = [arg.id '_' num2str(i)];
        target(i).location_id     = ['Location_' target(i).id ];
        target(i).orientation_id  = ['Orientation_' target(i).id ];
        target(i).scale_id        = ['Scale_' target(i).id ];
        target(i).model_id        = ['Model_' target(i).id ];
        
        placemark   = this.xml.createElement('Placemark');
        model       = this.xml.createElement('Model');
        location    = this.xml.createElement('Location');
        orientation = this.xml.createElement('Orientation');
        scale       = this.xml.createElement('Scale');
        link        = this.xml.createElement('Link');

        placemark.setAttribute('id',target(i).id);
        placemark.appendChild(this.textNode('name',arg.name));
        placemark.appendChild(this.textNode('visibility',num2str(arg.visibility)));
        placemark.appendChild(this.textNode('description',arg.description));

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
        
        location.setAttribute('id',target(i).location_id);
        location.appendChild(this.textNode('latitude',  num2str(lat(i),16)));
        location.appendChild(this.textNode('longitude', num2str(long(i),16)));
        location.appendChild(this.textNode('altitude',  num2str(alt(i),16)));

        orientation.setAttribute('id',target(i).orientation_id);
        orientation.appendChild(this.textNode('heading', num2str(heading(i))));
        orientation.appendChild(this.textNode('tilt',    num2str(tilt(i))));
        orientation.appendChild(this.textNode('roll',    num2str(roll(i))));   

        scale.setAttribute('id',target(i).scale_id);
        scale.appendChild(this.textNode('x', num2str(arg.scaleX(i))));
        scale.appendChild(this.textNode('y', num2str(arg.scaleY(i))));
        scale.appendChild(this.textNode('z', num2str(arg.scaleZ(i))));       

        link.appendChild(this.textNode('href',arg.model));

        model.setAttribute('id',target(i).model_id);
        model.appendChild(this.textNode('altitudeMode',arg.altitudeMode));

        model.appendChild(location);
        model.appendChild(orientation);
        model.appendChild(scale);
        model.appendChild(link);

        placemark.appendChild(model);
        this.doc.appendChild(placemark);
    end
   this.addIncludeFile(arg.model); 
    
end