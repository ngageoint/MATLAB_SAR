function target = quiver(this,long,lat,u,v,varargin)
%KML.QUIVER(long,lat,u,v) Create a quiver plot, similar to built-in
%   function quiver3, using line arrows. 
%  
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','','coordinates_type','','coordinates_id','');
    
    [longDEG,latDEG] = this.checkUnit(long,lat);
    
    p = inputParser;
    
    nlat = numel(lat);

    p.addRequired('lat',  @(a)isnumeric(a) && ~isempty(a));
    p.addRequired('long', @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('u',    @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('v',    @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);

    p.addParamValue('altitude',1000,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('scale',1,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('arrowBaseSize',0.3,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('plotArrows',true,@islogical);
    p.addParamValue('arrowHeadSize',0.3,@(a)isnumeric(a) && numel(a)==1);
    
    p.addParamValue('color','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
    p.addParamValue('name','kml_quiver',@ischar);
    p.addParamValue('id',kml.getTempID('kml_quiver'),@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,u,v,varargin{:});
    
    arg = p.Results;

    f = this.createFolder(arg.name);

    uv = sqrt((u+v).^2);
        
    edgelen = @(a)(max(a(:))-min(a(:))).^2;
    dS = edgelen(latDEG) + edgelen(longDEG);
    
    scale = sqrt((u.^2 + v.^2)./dS);
    scale = arg.scale .* scale;
    alpha = arg.arrowHeadSize;
    beta  = arg.arrowBaseSize;
    for i = 1:numel(lat)
        if arg.plotArrows
            long2 = scale(i).*[0 u(i) u(i)-alpha*(u(i)+beta*(v(i)+eps)) u(i) u(i)-alpha*(u(i)-beta*(v(i)+eps))] + long(i);
            lat2  = scale(i).*[0 v(i) v(i)-alpha*(v(i)-beta*(u(i)+eps)) v(i) v(i)-alpha*(v(i)+beta*(u(i)+eps))] + lat(i);
        else
            long2 = scale(i).*[0 u(i)] + long(i);
            lat2  = scale(i).*[0 v(i)] + lat(i);
        end
        target(i) = f.plot(long2,lat2, 'altitude',arg.altitude,...
                           'altitudeMode',arg.altitudeMode, ...
                           'visibility',arg.visibility, ...
                           'name',sprintf('Arrow %i',i), ...
                           'lineColor',arg.color, ...
                           'timeStamp', arg.timeStamp , ...
                           'timeSpanBegin', arg.timeSpanBegin , ...
                           'timeSpanEnd', arg.timeSpanEnd, ...                                                                  
                           'id', [arg.id '_' num2str(i)] ...
                           );
    end
end   