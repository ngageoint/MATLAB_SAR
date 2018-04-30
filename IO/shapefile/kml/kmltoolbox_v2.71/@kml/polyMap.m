function target = polyMap(this,long,lat,value,varargin)
%KML.POLY3(long,lat,alt) Draw a polygonal 3D form, with vertices given by long, lat and alt.
%  
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','');
    
    lat   = lat(:);
    long  = long(:);
    value = value(:);
    
    p = inputParser;
    nlat = numel(lat);

    p.addRequired('lat', @(a)isnumeric(a) &&~isempty(a));
    p.addRequired('long',@(a)isnumeric(a) &&~isempty(a) && numel(a)==nlat);
    p.addRequired('value', @(a)isnumeric(a) &&~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('id',kml.getTempID('kml_polyMap'),@ischar);
    p.addParamValue('name','kml_polyMap',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('colorMap','jet',@ischar);
    p.addParamValue('scale',1e4,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('width',1,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('type','square',@(a)ismember(lower(a),{'square','circle'}));
    p.addParamValue('noFolder',false,@islogical)
    
    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,value,varargin{:});
    
    arg = p.Results;
    
    fDR = this.checkUnit(1);

    lls = max(max(lat)-min(lat),max(long)-min(long));
    width = 0.5*lls./sqrt(numel(value));
    
    
    if arg.noFolder
        f = this;
    else
        f = this.createFolder(arg.name);
    end
        
    maxVal = max(value);
    minVal = min(value);
    
    if maxVal==minVal
        minVal = 0;
    end

    
    
    
    scale = arg.scale .* value; %-minVal)./(maxVal-minVal);
    
    ncolors = 100;
    cmap = feval(arg.colorMap,ncolors);
    vs = linspace(minVal,maxVal,ncolors);    

    switch lower(arg.type)
        case 'square'
            pLong = 0.5*[-1 -1 1 1 -1].*arg.width.*width;
            pLat  = 0.5*[-1 1 1 -1 -1].*arg.width.*width;
        case 'circle'
            t = linspace(0,2*pi,100);
            pLong = 0.5*sin(t).*arg.width.*width;
            pLat  = 0.5*cos(t).*arg.width.*width;
    end
    
    maxValCount = 1e4;
    
    
    if numel(value) > maxValCount
        warning('kml:maxPoly','Ignoring any point under the top %i values',maxValCount);
        [tmp,indexToPlot] = sort(-value);
        indexToPlot = indexToPlot(1:maxValCount).';
    else
        indexToPlot = 1:numel(value);
    end
    
    for i = indexToPlot
        if scale(i)>0
            color = [interp1(vs,cmap(:,1),value(i)) interp1(vs,cmap(:,2),value(i)) interp1(vs,cmap(:,3),value(i))];
            colorHex = kml.color2kmlHex(color);        
            target(end+1) = f.poly3(long(i) + pLong,lat(i) + pLat.*cosd(fDR*lat(i)),ones(size(pLong)).*scale(i), ...
                                    'extrude', true, ...
                                    'visibility',arg.visibility,...
                                    'altitudeMode','absolute',...
                                    'lineWidth',0, ...
                                    'polyColor',colorHex, ...
                                    'id', [arg.id '_' num2str(i)]);
        end
    end
    target(1) = []; %remove the empty initial field
end