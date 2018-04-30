function target = scatter(this,long,lat,varargin)
%KML.SCATTER(long,lat,alt) Places point markers in the positions given by long and lat
%  Similar to built-in scatter function. For a list of available markers, run
%    help kml.parseIconURL
%  
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    target = struct('type','','id','');

    [long,lat] = this.checkUnit(long,lat);

    nlat = numel(lat);

    p = inputParser;
    
    p.addRequired('lat', @(a)isnumeric(a) && ~isempty(a));
    p.addRequired('long',@(a)isnumeric(a) && ~isempty(a) && numel(a)==nlat);
    
    p.addParamValue('id',kml.getTempID('kml_scatter'),@ischar);
    p.addParamValue('name','kml_scatter',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('iconColor','FFFFFFFF',@(a)(iscell(a) && numel(a)==nlat) || (size(a,1)==1 || size(a,1)==nlat) || (ischar(a) && numel(a)==8));
    p.addParamValue('iconURL','',@(a)ischar(a)||(iscell(a)&&numel(a)==nlat));
    p.addParamValue('iconScale',1,@(a)isnumeric(a) && (numel(a)==1 || numel(a)==nlat));
    p.addParamValue('altitudeMode','clampToGround',@(a)ismember(a,{'clampToGround','relativeToGround','absolute'}));
    p.addParamValue('altitude',1,@(a)isnumeric(a) && (numel(a)==1 || numel(a)==nlat));
         
    p.addParamValue('timeStamp','',@(a)ischar(a)||iscell(a));
    p.addParamValue('timeSpanBegin','',@(a)ischar(a)||iscell(a));
    p.addParamValue('timeSpanEnd','',@(a)ischar(a)||iscell(a));
    
    p.parse(lat,long,varargin{:});
    
    arg = p.Results;

    if numel(arg.altitude)~=nlat
        alt = repmat(arg.altitude(1),size(lat));
    else
        alt = arg.altitude;
    end

    nc = numel(arg.iconColor);
    ns = numel(arg.iconScale);
    

    if ~iscell(arg.timeStamp)
        arg.timeStamp = {arg.timeStamp};
    end
    if numel(arg.timeStamp)~=nlat
        arg.timeStamp = repmat(arg.timeStamp,nlat,1);
    end
    
    if ~iscell(arg.timeSpanBegin)
        arg.timeSpanBegin = {arg.timeSpanBegin};
    end
    if numel(arg.timeSpanBegin)~=nlat
        arg.timeSpanBegin = repmat(arg.timeSpanBegin,nlat,1);
    end            

    if ~iscell(arg.timeSpanEnd)
        arg.timeSpanEnd = {arg.timeSpanEnd};
    end
    if numel(arg.timeSpanEnd)~=nlat
        arg.timeSpanEnd = repmat(arg.timeSpanEnd,nlat,1);
    end            
    
    if nc~=1
        if ~iscell(arg.iconColor) && ~ischar(arg.iconColor)
            if ismember(size(arg.iconColor,2),[3 4])
                %Color matrix, convert to string format
                iC = cell(size(arg.iconColor,1),1);
                for i = 1:size(arg.iconColor,1)
                    c = min(max(floor(arg.iconColor(i,:)*255),0),255);
                    if numel(c)==3
                        [r,g,b,a] = deal(c(1),c(2),c(3),255); 
                    else
                        [r,g,b,a] = deal(c(1),c(2),c(3),c(4)); 
                    end
                    [rhex, ghex, bhex, ahex ]= deal(dec2hex(r),dec2hex(g),dec2hex(b),dec2hex(a));
                    if length(rhex)==1,rhex=['0' rhex];end
                    if length(ghex)==1,ghex=['0' ghex];end
                    if length(bhex)==1,bhex=['0' bhex];end
                    if length(ahex)==1,ahex=['0' ahex];end
                    
                    iC{i} = [ahex bhex ghex rhex];
                end
                arg.iconColor = iC;
                nc = numel(arg.iconColor);
            else
                error('Invalid iconColor argument')
            end
        end
    end

    iconURL = kml.parseIconURL(arg.iconURL);
    scatterfolder = this.xml.createElement('Folder');    
    
    for i = 1:nlat
        target(i).id   = [arg.id '_' num2str(i)]; 
        target(i).type = 'Placemark';
        %         coordinates = mat2str([long(i) lat(i) alt(i)]);
        %         coordinates = coordinates(2:end-1);
        %         coordinates = strrep(coordinates,' ',',');
        %         coordinates = strrep(coordinates,';',' ');
        coordinates  = sprintf('%0.16g,%0.16g,%0.16g ',[long(i) lat(i) alt(i)].');    

        placemark   = this.xml.createElement('Placemark');
        point       = this.xml.createElement('Point');
        style       = this.xml.createElement('Style');
        iconstyle   = this.xml.createElement('IconStyle');
        icon        = this.xml.createElement('Icon');

        placemark.setAttribute('id',target(i).id);
        placemark.appendChild(this.textNode('name','')); %[arg.name '_' num2str(i)]
        placemark.appendChild(this.textNode('visibility',num2str(arg.visibility)));
        placemark.appendChild(this.textNode('description',arg.description));

        if ~isempty(arg.timeStamp{i})
            timeStamp = this.xml.createElement('TimeStamp');
            timeStamp.appendChild(this.textNode('when',arg.timeStamp{i}));
            placemark.appendChild(timeStamp);
        elseif ~isempty(arg.timeSpanBegin{i}) || ~isempty(arg.timeSpanEnd{i})
            timeSpan = this.xml.createElement('TimeSpan');
            if ~isempty(arg.timeSpanBegin{i})
                timeSpan.appendChild(this.textNode('begin',arg.timeSpanBegin{i}));
            end
            
            if ~isempty(arg.timeSpanEnd{i})
                timeSpan.appendChild(this.textNode('end',arg.timeSpanEnd{i}));
            end
            placemark.appendChild(timeSpan);
        end
        
        if nc==nlat
            iconstyle.appendChild(this.textNode('color',arg.iconColor{i}));
        else
            iconstyle.appendChild(this.textNode('color',arg.iconColor));
        end
        
        if ns==nlat
            iconstyle.appendChild(this.textNode('scale',num2str(arg.iconScale(i))));
        else
            iconstyle.appendChild(this.textNode('scale',num2str(arg.iconScale)));
        end
        
        if iscell(iconURL)
            icon.appendChild(this.textNode('href',iconURL{i}));
        else
            icon.appendChild(this.textNode('href',iconURL));
        end

        point.setAttribute('id',['Point_' arg.id]);
        point.appendChild(this.textNode('altitudeMode',arg.altitudeMode));
        point.appendChild(this.textNode('coordinates',coordinates));

        iconstyle.appendChild(icon);
        style.appendChild(iconstyle);
        placemark.appendChild(style);
        placemark.appendChild(point);
        scatterfolder.appendChild(placemark);
    end

    scatterfolder.appendChild(this.textNode('name',arg.name));
    scatterfolder.appendChild(this.textNode('id',  arg.id));  
    this.doc.appendChild(scatterfolder);
end