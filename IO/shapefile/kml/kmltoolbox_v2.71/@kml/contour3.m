function target = contour3(this,long,lat,alt,varargin)
%KML.CONTOUR(long,lat,alt) Create a 3D contour of alt in a grid defined by long and lat.
%  The altitude of each line is given by its associated contour level
%  Similar to built-in contour3 function
%
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','','coordinates_type','','coordinates_id','');
    
    p = inputParser;
    
    nlat = numel(lat);

    p.addRequired('lat',  @(a)isnumeric(a) && ~isempty(a));
    p.addRequired('long', @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('alt',  @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('name','kml_contour3',@ischar);
    p.addParamValue('id',kml.getTempID('kml_contour3'),@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('extrude',false,@islogical);
    p.addParamValue('colorMap','jet',@ischar);
    p.addParamValue('numberOfLevels','auto',@(a)(ischar(a) && strcmpi(a,'auto')) ||(isnumeric(a)));
    p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));

    p.addParamValue('noFolder',false,@islogical)
    
    p.addParamValue('showText',false,@islogical)
    p.addParamValue('levelStep',1,@isnumeric)
    p.addParamValue('labelSpacing',inf,@isnumeric)

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    

    p.addParamValue('lineWidth',1,@isnumeric);
    
    p.parse(lat,long,alt,varargin{:});
    
    arg = p.Results;
    
    if arg.noFolder
        f = this;
    else
        f = this.createFolder(arg.name);
    end
    
    if isnumeric(arg.numberOfLevels)
        c = contours(lat,long,alt,arg.numberOfLevels);
    else
        c = contours(lat,long,alt);
    end
    
    minAlt = min(min(alt));
    maxAlt = max(max(alt));
    
    i = 1; k=1;
    while true
        l = c(1,i);
        n = c(2,i);
        C(k).lat = c(1,(i+1):(i+n));
        C(k).long = c(2,(i+1):(i+n));
        C(k).l = l;
        i = i + n + 1;
        k = k + 1;
        if i > size(c,2)
            break;
        end
    end
    
    levels = unique([C(:).l]);
    shownLevels = levels(1:arg.levelStep:numel(levels));
    
    ncolors = 100;
    cmap = feval(arg.colorMap,ncolors);
    levelsCMAP = linspace(min(levels),max(levels),ncolors);
    middleLevel = levels(ceil(numel(levels)/2));
    
    for i = 1:numel(C)
        
        clev = C(i).l;
        iC = find(levels==clev);
        if clev <= middleLevel
            clev = clev - 1;
        end
        
        if numel(levels) > 1
            iC = round(interp1(levelsCMAP,linspace(0,ncolors-1,numel(levelsCMAP)),clev,'linear',0));
        else
            iC = 1;
        end        
        
        color = cmap(iC+1 ,:);

        colorHex = kml.color2kmlHex(color);

        polyColor = colorHex;
        
        if arg.extrude
            polyColor(1:2) = 'FF';
            colorHex(1:2) = '00';
            target(end+1) = f.poly(C(i).long,C(i).lat, 'lineColor', colorHex, ...
                                       'altitudeMode',arg.altitudeMode, ...
                                       'altitude',C(i).l,...
                                       'visibility',arg.visibility, ...
                                       'name',sprintf('Level %g',C(i).l), ...
                                       'timeStamp', arg.timeStamp , ...
                                       'timeSpanBegin', arg.timeSpanBegin , ...
                                       'timeSpanEnd', arg.timeSpanEnd, ...     
                                       'extrude',arg.extrude, ...
                                       'polyColor',polyColor, ...
                                       'id',[arg.id '_poly_' num2str(i)] ...  
                                       );
        else
            polyColor(1:2) = '00';
            
            target(end+1) = f.plot(C(i).long,C(i).lat, 'lineColor', colorHex, ...
                                      'lineWidth',arg.lineWidth,...
                                       'altitudeMode',arg.altitudeMode, ...
                                       'altitude',C(i).l,...
                                       'visibility',arg.visibility, ...
                                       'name',sprintf('Level %g',C(i).l), ...
                                       'timeStamp', arg.timeStamp , ...
                                       'timeSpanBegin', arg.timeSpanBegin , ...
                                       'timeSpanEnd', arg.timeSpanEnd, ...     
                                       'extrude',arg.extrude, ...
                                       'polyColor',polyColor, ...
                                       'id',[arg.id '_plot_' num2str(i)] ...  
                                       );
        end
        
        if arg.showText && ismember(C(i).l,shownLevels)
            N = numel(C(i).lat);
            altTxt = C(i).l;
            if ~isfinite(arg.labelSpacing)
                R = randi(N,1);
                latTxt  = interp1(1:N,C(i).lat,R);
                longTxt = interp1(1:N,C(i).long,R);
                
                target(end+1) = f.text(longTxt,latTxt,altTxt,sprintf('%g',C(i).l), 'id',[arg.id '_text_' num2str(i)]);
            else
               dist = kml.ll2dist(C(i).long(1:end-1),C(i).lat(1:end-1), ...
                                  C(i).long(2:end),C(i).lat(2:end));
               dist = [0 cumsum(dist)];
               for d = 1:arg.labelSpacing:dist(end);
                   latTxt  = interp1(dist,C(i).lat,d);
                   longTxt = interp1(dist,C(i).long,d);
                   target(end+1) = f.text(longTxt,latTxt,altTxt,sprintf('%g',C(i).l), 'id',[arg.id '_text_' num2str(i) '_' num2str(d)]);
               end
            end
        end
        
    end
    target(1) = []; %remove the empty initial field
end