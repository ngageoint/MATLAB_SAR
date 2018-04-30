function target = surf(this,long,lat,alt,varargin)
%KML.surf(long,lat,alt) Create a surface contour of alt in a grid defined by long and lat. 
%   Similar to built-in surf function
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','');

    p = inputParser;
    
    [long,lat] = this.checkUnit(long,lat);

    
    nlat = numel(lat);
    [r,c] = size(lat);
    p.addRequired('lat',  @(a)isnumeric(a) && ~isempty(a));
    p.addRequired('long', @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('alt',  @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('name','kml_surf',@ischar);
    p.addParamValue('id',kml.getTempID('kml_surf'),@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('colorMap','jet',@ischar);
    p.addParamValue('color','',@(a)ischar(a) && numel(a)==8);

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,alt,varargin{:});
    
    arg = p.Results;
    arg.visibility = num2str(arg.visibility);
    
    if isempty(arg.color)
        lineColor = 'FF000000';
    else
        lineColor = '00000000'; %arg.color;
    end
    
    
    f = this.createFolder(arg.name);

    ncolors = 100;
    cmap = feval(arg.colorMap,ncolors);
    cspace = 1:ncolors;
    aspace = linspace(min(alt(:)),max(alt(:)),ncolors);
    for i = 2:r
        for j = 2:c
           longC = [long(i-1,j-1) long(i,j-1) long(i,j) long(i-1,j-1)];% long(i-1,j) long(i,j) long(i-1,j-1)];   
           latC  = [lat(i-1,j-1)  lat(i,j-1)  lat(i,j)  lat(i-1,j-1)];%  lat(i-1,j) lat(i,j) lat(i-1,j-1)];
           altC  = [alt(i-1,j-1)  alt(i,j-1)  alt(i,j)  alt(i-1,j-1)];%  alt(i-1,j) alt(i,j) alt(i-1,j-1)];
           
           alev = mean([alt(i-1,j-1) alt(i,j-1)  alt(i,j)]);
           iC = round(interp1(aspace,cspace,alev,'linear',1));
           color = cmap(iC ,:);
           
           if isempty(arg.color)
               colorHex = kml.color2kmlHex(color);
           else
               colorHex = arg.color;
           end
           
           
%             f.poly3(longC,latC,altC, 'polyColor', colorHex, ...
%                                        'lineColor','FF000000',... %['00' colorHex(3:end)],...
%                                        'altitudeMode','relativeToGround', ...
%                                        'visibility',arg.visibility, ...
%                                        'name',sprintf('Cell (%i,%i)',i,j), ...
%                                        'timeStamp', arg.timeStamp , ...
%                                        'timeSpanBegin', arg.timeSpanBegin , ...
%                                        'timeSpanEnd', arg.timeSpanEnd ...
%                                        );

           target(end+1).id = fastPoly(longC,latC,altC,sprintf('Cell (%i,%i)',i,j),[arg.id '_' sprintf('Cell (%i,%i)',i,j)]);
           target(end).type = 'Placemark';
    
           longC = [long(i-1,j-1) long(i-1,j) long(i,j) long(i-1,j-1)];
           latC  = [lat(i-1,j-1)  lat(i-1,j)  lat(i,j)  lat(i-1,j-1)];
           altC  = [alt(i-1,j-1)  alt(i-1,j)  alt(i,j)  alt(i-1,j-1)];
           
           alev  = mean([alt(i-1,j-1)  alt(i-1,j)  alt(i,j)]);
           iC    = round(interp1(aspace,cspace,alev,'linear',1));
           color = cmap(iC ,:);
           
           if isempty(arg.color)
               colorHex = kml.color2kmlHex(color);
           else
               colorHex = arg.color;
           end
           
           
%             f.poly3(longC,latC,altC, 'polyColor', colorHex, ...
%                                        'lineColor',lineColor,... %['00' colorHex(3:end)],...
%                                        'altitudeMode','relativeToGround', ...
%                                        'visibility',arg.visibility, ...
%                                        'name',sprintf('Cell (%i,%i)',i,j), ...
%                                        'timeStamp', arg.timeStamp , ...
%                                        'timeSpanBegin', arg.timeSpanBegin , ...
%                                        'timeSpanEnd', arg.timeSpanEnd ...
%                                        );
% 

            target(end+1).id = fastPoly(longC,latC,altC,sprintf('Cell 2 (%i,%i)',i,j),[arg.id '_' sprintf('Cell (%i,%i)',i,j)]);
            target(end).type = 'Placemark';
        end
    end
    
    
    function id = fastPoly(long,lat,alt,name,id)

        extrudeNode      = this.textNode('extrude','0');
        tessellateNode   = this.textNode('tesselate','1');
        altitudeModeNode = this.textNode('altitudeMode','absolute');
        visibilityNode   = this.textNode('visibility',arg.visibility);
        widthNode        = this.textNode('width','1');
        lineColorNode    = this.textNode('color',lineColor);
        
        coordinates  = sprintf('%0.16g,%0.16g,%0.16g ',[long(:) lat(:) alt(:)].');    

        placemark   = this.xml.createElement('Placemark');

        polygon     = this.xml.createElement('Polygon');
        outboundary = this.xml.createElement('outerBoundaryIs');
        linearring  = this.xml.createElement('LinearRing');
        style       = this.xml.createElement('Style');
        linestyle   = this.xml.createElement('LineStyle');
        polystyle   = this.xml.createElement('PolyStyle');

        placemark.setAttribute('id',id);
        placemark.appendChild(this.textNode('name',name));
        placemark.appendChild(visibilityNode);
    
        linestyle.appendChild(lineColorNode);
        linestyle.appendChild(widthNode);

        polystyle.appendChild(this.textNode('color',colorHex));

        linearring.setAttribute('id','LinearRing');
        linearring.appendChild(this.textNode('coordinates',coordinates));

        polygon.setAttribute('id','Polygon');
        polygon.appendChild(extrudeNode);
        polygon.appendChild(tessellateNode);
        polygon.appendChild(altitudeModeNode);

        outboundary.appendChild(linearring);
        polygon.appendChild(outboundary);

        style.appendChild(linestyle);
        style.appendChild(polystyle);
        placemark.appendChild(style);
        placemark.appendChild(polygon);
        f.doc.appendChild(placemark);
    end
end