function target = quiver3d(this,long,lat,alt,u,v,w,varargin)
%KML.QUIVER3D(long,lat,alt,u,v,w) Create a quiver plot, similar to built-in
%   function quiver3, using 3D arrows. 
%  
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','','location_id','','orientation_id','','scale_id','','model_id','');
    
    [longDEG,latDEG] = this.checkUnit(long,lat);
    
    p = inputParser;
    
    nlat = numel(lat);

    p.addRequired('lat',  @(a)isnumeric(a) && ~isempty(a));
    p.addRequired('long', @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('alt',  @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('u',    @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('v',    @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    p.addRequired('w',    @(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);

    p.addParamValue('scale',100,@(a)isnumeric(a) && numel(a)==1);
    p.addParamValue('color','FFFFFFFF',@(a)ischar(a) && numel(a)==8);
    p.addParamValue('name','kml_quiver3D',@ischar);
    p.addParamValue('id',kml.getTempID('kml_quiver3D'),@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('model','',@ischar);
    p.addParamValue('altitudeMode','relativeToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    p.parse(lat,long,alt,u,v,w,varargin{:});
    
    arg = p.Results;

    if isempty(arg.model)
        arrowfile = 'arrow3d.dae';

        arg.model = arrowfile;
    end
    
    modelpath = which(arg.model);

    if ~isempty(modelpath)
        if isempty(dir(arg.model))
            copyfile(modelpath,arg.model);
        end
    else
        error('File %s not found!',arg.model);
    end
        
    f = this.createFolder(arg.name);

    uv = sqrt((u+v).^2);
    if strcmpi(this.unit,'deg')
        rad2deg = 180/pi;
    else
        rad2deg = 1;
    end
    heading = rad2deg * atan2(u,v);
    tilt    = rad2deg * atan2(w,uv)-pi/2;
    roll    = zeros(size(lat));
    
    
    edgelen = @(a)(max(a(:))-min(a(:))).^2;
    dS = edgelen(latDEG) + edgelen(longDEG) + edgelen(alt);
    
    scale = sqrt((u.^2 + v.^2 + w.^2)./dS);
    scale = arg.scale .* scale;
    
    for i = 1:numel(lat)
        target(i) = f.model(long(i),lat(i),alt(i),heading(i),tilt(i),roll(i), 'scale',scale(i),'model',arg.model, ...
                                                                  'altitudeMode',arg.altitudeMode, ...
                                                                  'visibility',arg.visibility, ...
                                                                  'name',sprintf('Arrow %i',i), ...
                                                                  'timeStamp', arg.timeStamp , ...
                                                                  'timeSpanBegin', arg.timeSpanBegin , ...
                                                                  'timeSpanEnd', arg.timeSpanEnd, ...        
                                                                  'id', [arg.id '_' num2str(i)] ...
                                                                  );
    end
    
    this.addIncludeFile(arg.model);
end   