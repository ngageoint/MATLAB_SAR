function target = quadoverlay(this,long,lat,varargin)
%KML.quadoverlay(longitude, latitude) Places the image file (specified by the pair 
%  attribute 'file','image.png') as a ground overlay in the kml file. The corners  
%  of the image are given by inputs longitude and latitude (which should represent 4 points). 
%  To make the overlay transparent, change the alpha portion of the color parameter
%  to a different hex value - eg.: 50% transparent, use KML.OVERLAY(...,'color','80FFFFFF')
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','');

    [long,lat] = this.checkUnit(long,lat);
    
    p = inputParser;
       
    p.addRequired('long', @(a)isnumeric(a) && numel(a)==4);
    p.addRequired('lat', @(a)isnumeric(a) && numel(a)==4);
    
    p.addParamValue('id',kml.getTempID('kml_overlay'),@ischar);
    p.addParamValue('name','kml_overlay',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('file','',@ischar);
    p.addParamValue('viewBoundScale',1,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('color','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
    p.addParamValue('altitude',1,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);
    p.addParamValue('drawOrder',0,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);
    p.addParamValue('altitudeMode','clampToGround',@(a)ismember(a,{'clampToGround','absolute'}));

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(long,lat,varargin{:});
    
    arg = p.Results;

    if isempty(arg.file)
        error('Missing file parameter')
    end
    
    coordinates  = sprintf('%0.16g,%0.16g ',[long(:) lat(:)].');    
    
    overlay      = this.xml.createElement('GroundOverlay');
    icon         = this.xml.createElement('Icon');
    gxLatLonQuad = this.xml.createElement('gx:LatLonQuad');
    
    overlay.setAttribute('id',arg.id);
    overlay.appendChild(this.textNode('name',arg.name));
    overlay.appendChild(this.textNode('color',arg.color));
    overlay.appendChild(this.textNode('visibility',num2str(arg.visibility)));
    overlay.appendChild(this.textNode('description',arg.description));
    overlay.appendChild(this.textNode('altitude',num2str(arg.altitude)));
    overlay.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
    overlay.appendChild(this.textNode('drawOrder',num2str(arg.drawOrder)));

    if ~isempty(arg.timeStamp)
        timeStamp = this.xml.createElement('TimeStamp');
        timeStamp.appendChild(this.textNode('when',arg.timeStamp));
        overlay.appendChild(timeStamp);
    elseif ~isempty(arg.timeSpanBegin) || ~isempty(arg.timeSpanEnd)
        timeSpan = this.xml.createElement('TimeSpan');
        if ~isempty(arg.timeSpanBegin)
            timeSpan.appendChild(this.textNode('begin',arg.timeSpanBegin));
        end
        
        if ~isempty(arg.timeSpanEnd)
            timeSpan.appendChild(this.textNode('end',arg.timeSpanEnd));
        end
        overlay.appendChild(timeSpan);
    end    
    
    gxLatLonQuad.appendChild(this.textNode('coordinates',coordinates));

    icon.appendChild(this.textNode('href',arg.file));
    icon.appendChild(this.textNode('viewBoundScale',num2str(arg.viewBoundScale)));

    overlay.appendChild(gxLatLonQuad);
    overlay.appendChild(icon);
    this.doc.appendChild(overlay);
    
    target.id   = arg.id;
    target.type = 'GroundOverlay';
    
    this.addIncludeFile(arg.file);
end