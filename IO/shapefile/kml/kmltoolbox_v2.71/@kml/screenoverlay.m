function target = screenoverlay(this,overlayXY,screenXY,sz,varargin)
%KML.SCREENOVERLAY(overlayXY,screenXY,sz) Places the image file (specified by the pair 
%  attribute 'file','image.png') as a screen overlay in the kml file.
%  Check https://developers.google.com/kml/documentation/kmlreference#screenoverlay for more info.
%  To make the overlay transparent, change the alpha portion of the color parameter
%  to a different hex value - eg.: 50% transparent, use KML.OVERLAY(...,'color','80FFFFFF')
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.4 $  $Date: 2012/12/07 08:00:00 $

    target = struct('type','','id','');

  
    
    p = inputParser;
       
    p.addRequired('overlayXY', @(a)isnumeric(a) && ~isempty(a) && numel(a)==2);
    p.addRequired('screenXY', @(a)isnumeric(a) && ~isempty(a) && numel(a)==2);
    p.addRequired('sz',@(a)isnumeric(a) && ~isempty(a) && numel(a)==2);
    
    p.addParamValue('id',kml.getTempID('kml_overlay'),@ischar);
    p.addParamValue('name','kml_overlay',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('file','',@ischar);
    p.addParamValue('viewBoundScale',1,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('color','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
    p.addParamValue('altitude',1,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);
    p.addParamValue('rotation',0,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);
    p.addParamValue('drawOrder',0,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);
    p.addParamValue('altitudeMode','clampToGround',@(a)ismember(a,{'clampToGround','absolute'}));
    p.addParamValue('units','fraction',@(a)ismember(a,{'fraction','pixels','insetPixels'}));

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(overlayXY,screenXY,sz,varargin{:});
    
    arg = p.Results;

    if isempty(arg.file)
        error('Missing file parameter')
    end
    
    arg.rotation = this.checkUnit(arg.rotation);
    
    arg.overlayPos = overlayXY;
    arg.screenPos  = screenXY;
    arg.sz         = sz;
    
% <ScreenOverlay id="khScreenOverlay756">
%   <name>Simple crosshairs</name>
%   <description>This screen overlay uses fractional positioning
%    to put the image in the exact center of the screen</description>
%   <Icon>
%     <href>http://myserver/myimage.jpg</href>
%   </Icon>
%   <overlayXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
%   <screenXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
%   <rotation>39.37878630116985</rotation>
%   <size x="0" y="0" xunits="pixels" yunits="pixels"/>
% </ScreenOverlay>


%%

    overlay    = this.xml.createElement('ScreenOverlay');
    overlayXY  = this.xml.createElement('overlayXY');
    screenXY   = this.xml.createElement('screenXY');
    size       = this.xml.createElement('size');
    icon       = this.xml.createElement('Icon');

    overlay.setAttribute('id',arg.id);
    overlay.appendChild(this.textNode('name',arg.name));
    overlay.appendChild(this.textNode('open','1'));
    overlay.appendChild(this.textNode('visibility',num2str(arg.visibility)));
    overlay.appendChild(this.textNode('description',arg.description));
    overlay.appendChild(this.textNode('drawOrder',num2str(arg.drawOrder)));
    overlay.appendChild(this.textNode('color',arg.color));

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
    
    overlayXY.setAttribute('x',num2str(arg.overlayPos(1),16))
    overlayXY.setAttribute('y',num2str(arg.overlayPos(2),16))
    overlayXY.setAttribute('xunits','fraction')
    overlayXY.setAttribute('yunits','fraction')
    
    screenXY.setAttribute('x',num2str(arg.screenPos(1),16))
    screenXY.setAttribute('y',num2str(arg.screenPos(2),16))
    screenXY.setAttribute('xunits','fraction')
    screenXY.setAttribute('yunits','fraction')    
    
    size.setAttribute('x',num2str(arg.sz(1)))
    size.setAttribute('y',num2str(arg.sz(2)))
    size.setAttribute('xunits','fraction')
    size.setAttribute('yunits','fraction')        
    overlay.appendChild(this.textNode('rotation', num2str(arg.rotation)));

    icon.appendChild(this.textNode('href',arg.file));

    overlay.appendChild(overlayXY);
    overlay.appendChild(screenXY);
    overlay.appendChild(size);
    overlay.appendChild(icon);
    
    this.doc.appendChild(overlay);   
    
    target.id   = arg.id;
    target.type = 'ScreenOverlay';
    this.addIncludeFile(arg.file);
end