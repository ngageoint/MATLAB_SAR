function target = colorbar(this,cAxis,colorMap,varargin)
    ss = get(0,'ScreenSize')*.8;

    p = inputParser;
    
    p.addRequired('cAxis', @isnumeric);
    p.addRequired('colorMap', @isnumeric);

    p.addParamValue('type','vertical',@(a)ismember(lower(a),{'vertical','horizontal'}));
    p.addParamValue('tick',[],@isnumeric);
    p.addParamValue('tickLabel',{},@iscell);
        
    

    p.addParamValue('id','kml_colorbar',@ischar);
    p.addParamValue('name','kml_colorbar',@ischar);
    p.addParamValue('description','',@ischar);
    p.addParamValue('visibility',true,@islogical);
    p.addParamValue('rotation',0,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);
    p.addParamValue('drawOrder',0,@(a)isnumeric(a) &&~isempty(a) && numel(a)==1);

    p.addParamValue('timeStamp','',@ischar);
    p.addParamValue('timeSpanBegin','',@ischar);
    p.addParamValue('timeSpanEnd','',@ischar);    
    
    
    p.parse(cAxis,colorMap,varargin{:});
    
    arg = p.Results;

    if numel(cAxis) > 2
        cAxis = [min(cAxis) max(cAxis)];
    end
    
    l = linspace(cAxis(1),cAxis(2),1000);
    fh = figure;
    
    
    if strcmpi(arg.type,'vertical')
        ratio = 40;
        
        pcolor([0*l.' 1+0*l.'],[l.' l.'],[l.' l.'])
        
        set(gcf,'Position',[0 0 ss(4)./ratio  ss(4)]);%[ss(1)+ss(4)./ratio  ss(2) ss(4)./ratio  ss(4)]);
        set(gca,'XTick',[]);
        if ~isempty(arg.tick)
            set(gca,'YTick',arg.tick);
        end
        if ~isempty(arg.tickLabel)
            set(gca,'YTickLabel',arg.tickLabel);
        end        
        sz         = [2/ratio 1];
        overlayPos = [2 1];
        screenPos  = [1 1];
    else
        ratio = 40;
        pcolor([l.' l.'],[0*l.' 1+0*l.'],[l.' l.'])
        set(gcf,'Position',[0 0  ss(3) ss(3)./ratio ]);%[ss(1) ss(2)+ss(3)./ratio  ss(3) ss(3)./ratio ]);
        set(gca,'YTick',[]);
        if ~isempty(arg.tick)
            set(gca,'XTick',arg.tick);
        end
        if ~isempty(arg.tickLabel)
            set(gca,'XTickLabel',arg.tickLabel);
        end                
        sz         = [1 2/ratio];
        overlayPos = [1 (ratio-4)/2];
        screenPos  = [1 1];        
    end
    
    movegui(fh,'north');
    
    shading flat;
    caxis(cAxis);
    colormap(colorMap);
    
    
    
    name = 'Colorbar';
    basename = 'Colorbar';
    
    name = [name '.png'];
    k = 1;
    while ~isempty(dir(name))
        name = sprintf('%s (%i)%s',basename,k,'.png');
        k = k+1;
    end
    
    bgColor = [0 0 0];
    fgColor = [1 1 1];
    
    set(gca,'Color',bgColor,'XColor',fgColor,'YColor',fgColor);
    set(fh,'Color',bgColor);
    set(gca,'FontWeight','bold','FontSize',15)
%     grid on
%     set(gca,'GridLineStyle','-')
    
%     figure(fh)
%     drawnow
%     pause(1) %drawing on windows starts with a transparent window, better give some time for it to show, so we don't capture the background
%     drawnow
%     drawnow('expose');
%     im = getframe(fh); 
%     drawnow('expose');
%     im = im.cdata;
    
    set(fh,'PaperPositionMode','auto','InvertHardcopy','off')
    print('-dpng','-r0',name,fh)
    im = imread(name);

    axc = ceil(bgColor*255);
    alphaMap = uint8(~(im(:,:,1)==axc(1) & im(:,:,2)==axc(2) & im(:,:,3)==axc(3))*255);

    imwrite(im,name,'Alpha',alphaMap);
    close(fh);
    drawnow


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
    
    overlayXY.setAttribute('x',num2str(overlayPos(1),16))
    overlayXY.setAttribute('y',num2str(overlayPos(2),16))
    overlayXY.setAttribute('xunits','fraction')
    overlayXY.setAttribute('yunits','fraction')
    
    screenXY.setAttribute('x',num2str(screenPos(1),16))
    screenXY.setAttribute('y',num2str(screenPos(2),16))
    screenXY.setAttribute('xunits','fraction')
    screenXY.setAttribute('yunits','fraction')    
    
    size.setAttribute('x',num2str(sz(1),16))
    size.setAttribute('y',num2str(sz(2),16))
    size.setAttribute('xunits','fraction')
    size.setAttribute('yunits','fraction')        
    overlay.appendChild(this.textNode('rotation', num2str(arg.rotation)));

    icon.appendChild(this.textNode('href',name));

    overlay.appendChild(overlayXY);
    overlay.appendChild(screenXY);
    overlay.appendChild(size);
    overlay.appendChild(icon);
    
    this.doc.appendChild(overlay);   
    
    target.id   = arg.id;
    target.type = 'ScreenOverlay';
    
    this.addIncludeFile(name);
end