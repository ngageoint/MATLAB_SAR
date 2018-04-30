function target = poly3(this,long,lat,alt,varargin)
%KML.POLY3(long,lat,alt) Draw a polygonal 3D form, with vertices given by long, lat and alt.
%  
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $
  
    target = struct('type','','id','');
  
    [long,lat] = this.checkUnit(long,lat);

    p = inputParser;
    nlat = numel(lat);

    p.addRequired('lat', @(a)isnumeric(a) && isvector(a) &&~isempty(a));
    p.addRequired('long',@(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('alt', @(a)isnumeric(a) && isvector(a) &&~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('id',kml.getTempID('kml_poly3'),@ischar);
    p.addParamValue('name','kml_poly3',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('lineColor','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
    p.addParamValue('polyColor','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
    p.addParamValue('lineWidth',1,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));
    p.addParamValue('extrude',false,@islogical);
    p.addParamValue('tessellate',true,@islogical);

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,alt,varargin{:});
    
    arg = p.Results;
    
%     coordinates = mat2str([long(:) lat(:) alt(:)]);
%     coordinates = coordinates(2:end-1);
%     coordinates = strrep(coordinates,' ',',');
%     coordinates = strrep(coordinates,';',' ');

    coordinates  = sprintf('%0.16g,%0.16g,%0.16g ',[long(:) lat(:) alt(:)].');    
    
    placemark   = this.xml.createElement('Placemark');
    
    polygon     = this.xml.createElement('Polygon');
    outboundary = this.xml.createElement('outerBoundaryIs');
    linearring  = this.xml.createElement('LinearRing');
    style       = this.xml.createElement('Style');
    linestyle   = this.xml.createElement('LineStyle');
    polystyle   = this.xml.createElement('PolyStyle');

    placemark.setAttribute('id',arg.id);
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
    
    linestyle.appendChild(this.textNode('color',arg.lineColor));
    linestyle.appendChild(this.textNode('width',num2str(arg.lineWidth)));
    
    polystyle.appendChild(this.textNode('color',arg.polyColor));
    
    linearring.setAttribute('id',['LinearRing_' arg.id]);
    linearring.appendChild(this.textNode('coordinates',coordinates));

    polygon.setAttribute('id',['Polygon_' arg.id]);
    polygon.appendChild(this.textNode('extrude',num2str(arg.extrude)));
    polygon.appendChild(this.textNode('tesselate',num2str(arg.tessellate)));
    polygon.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
    

    outboundary.appendChild(linearring);
    polygon.appendChild(outboundary);
    
    style.appendChild(linestyle);
    style.appendChild(polystyle);
    placemark.appendChild(style);
    placemark.appendChild(polygon);
    this.doc.appendChild(placemark);
    
    target.type = 'Placemark';
    target.id   = arg.id;
end