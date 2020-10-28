function ApToolnewpos(pos, hObject)
handles = guidata(hObject);

%get filter type
if get(handles.None,'Value')
    WindowFilter = 0;
elseif get(handles.Gaussian,'Value')
    WindowFilter = 1;
elseif get(handles.x4,'Value')
    WindowFilter = 2;
elseif get(handles.Hamming,'Value')
    WindowFilter = 3;
else
    WindowFilter = 4;
end
AlwaysApply = get(handles.AlwaysApply,'Value');

% Get min and max values
% Constraint function should avoid values outside bounds but we check anyway
xmin = max(handles.ZpLimsAz(1),floor(pos(1)));
xmax = min(handles.ZpLimsAz(2),floor(pos(1)+pos(3)-1));
ymin = max(handles.ZpLimsRn(1),floor(pos(2)));
ymax = min(handles.ZpLimsRn(2),floor(pos(2)+pos(4)-1));

%compute Az/Range resolution
CRFrac = ((xmax-xmin+1)/handles.ZpWidthAz)*100;
RangeFrac = ((ymax-ymin+1)/handles.ZpWidthRn)*100;

%compute resolutions
if get(handles.InverseSelection,'Value')
    RangeRes = handles.meta.Grid.Row.ImpRespWid./((100-RangeFrac)/100);
    AzRes = handles.meta.Grid.Col.ImpRespWid./((100-CRFrac)/100);
else
    RangeRes = handles.meta.Grid.Row.ImpRespWid./(RangeFrac/100);
    AzRes = handles.meta.Grid.Col.ImpRespWid./(CRFrac/100);
end

%set resolutions in the specified units
%azimuth
if get(handles.MetricCheck,'Value')
    %if less than a meter, display in cm
    if AzRes<1
        set(handles.AzRes,'String',round(AzRes*100));
        set(handles.AzUnits,'String','cm');
    else
        set(handles.AzRes,'String',roundn(AzRes,-1));
        set(handles.AzUnits,'String','meters');
    end
    if RangeRes<1
        set(handles.RgRes,'String',round(RangeRes*100));
        set(handles.RgUnits,'String','cm');
    else
        set(handles.RgRes,'String',roundn(RangeRes,-1));
        set(handles.RgUnits,'String','meters');
    end
else
    %convert to inches
    AzRes = 12*(AzRes./0.3048);
    RangeRes = 12*(RangeRes./0.3048);
    if AzRes<12
        set(handles.AzRes,'String',roundn(AzRes,-1));
        set(handles.AzUnits,'String','inches');
    else
        set(handles.AzRes,'String',roundn(AzRes/12,-1));
        set(handles.AzUnits,'String','feet');
    end
    if RangeRes<12
        set(handles.RgRes,'String',roundn(RangeRes,-1));
        set(handles.RgUnits,'String','inches');
    else
        set(handles.RgRes,'String',roundn(RangeRes/12,-1));
        set(handles.RgUnits,'String','feet');
    end
end

%compute time and bandwidth
[CollectTime,Bandwidth] = ComputeTimeBandwidth(handles.meta);
set(handles.ApTime,'String',roundn(CollectTime*CRFrac/100,-2));
set(handles.ApBandwidth,'String',round(Bandwidth*RangeFrac/1e8)); 

%need to modify xposition if flight is left
[ny, nx, nz] = size(handles.phasehistory);
if isfield(handles.meta,'SCPCOA') && isfield(handles.meta.SCPCOA,'SideOfTrack') && ...
        strcmp(handles.meta.SCPCOA.SideOfTrack,'L')
    xminold = xmin;
    xmin = round(nx-xmax);
    if (xmin<1);xmin=1;end
    xmax = round(nx-xminold);    
end

%now reform image with selected aperture
if get(handles.InverseSelection,'Value')
    if get(handles.InversePolar,'Value')
        %set inverse AOI base on k_a/k_r mapping to k_u/k_v space
        Ap = ComputeInvPolarApData([xmin xmax],[ymin ymax],handles);
    else
        Ap = handles.phasehistory;
        Ap(ymin:ymax,xmin:xmax,:) = 0;
    end
else
    if get(handles.InversePolar,'Value')
        %set inverse AOI base on k_a/k_r mapping to k_u/k_v space
        Ap = ComputeInvPolarApData([xmin xmax],[ymin ymax],handles);
    else
        Ap = zeros(ny,nx,nz);
        
        phdchip = handles.phasehistory(ymin:ymax,xmin:xmax,:);
           
        if (ymax-ymin) < handles.ZpWidthRn - 2
            FilterRange = 1;
        else
            FilterRange = 0;
        end
        if (xmax-xmin) < handles.ZpWidthAz - 2
            FilterAz = 1;
        else
            FilterAz = 0;
        end

        if AlwaysApply
            %apply filtering on both dimensions
            FilterRange = 1;
            FilterAz = 1;
        end
        
        
        phdchip = ApToolFilterPHDChip(phdchip,WindowFilter,FilterRange,FilterAz);

        Ap(ymin:ymax,xmin:xmax,:) = phdchip;    
    end
end

[ImMag,handles] = ApToolmakeDisplayable(handles,Ap);

remaps = get(handles.RemapCombo,'String');
remap = remaps{get(handles.RemapCombo,'Value')};
if strcmp(remap,'linearremap')
    % Rescale resolution to sigmas above mean from previous frame. This is
    % necessary since we lose energy if we reduce the aperture.
    chipmag = get(handles.imghandle,'CData');
    if isfloat(chipmag) % Previous frame was also linearremap
        CLim = get(handles.image,'CLim');
        old_sigma_over_mean = (CLim(2)-mean(chipmag(:)))./std(chipmag(:));
    else % Default value when coming from another remap type
        old_sigma_over_mean = 3;
    end
    CLim = [0 mean(single(ImMag(:)))+std(single(ImMag(:)))*old_sigma_over_mean];
    handles.apiSP1.replaceImage(ImMag,CLim,'PreserveView',true);
else
    handles.apiSP1.replaceImage(ImMag,'PreserveView',true);
end

if (handles.Record == 1)
    %record frame
    temp.Pos = pos;
    temp.CLim = get(handles.image,'CLim');
    temp.Mag = handles.apiSP1.getMagnification();
    R = handles.apiSP1.getVisibleImageRect();
    temp.Center = [round(R(1)+R(3)/2) round(R(2)+R(4)/2)]; 
    %save mode
    if isfield(handles, 'InLoop') && handles.InLoop
        %check for ResMode
        %if get(handles.ResModeCheck,'Value')
        %    temp.Mode = 'Multi-Resolution';
        %else
            %determine if this is slow of fast-time
            if get(handles.Slow,'Value')
                temp.Mode = 'Slow-Time';
            else
                temp.Mode = 'Fast-Time';
            end
        %end
    else
        temp.Mode = 'User Interactive';
    end
        
    %save resolution (Ground Plane)
    temp.AzRes = str2double(get(handles.AzRes,'String'));
    temp.RangeRes = str2double(get(handles.RgRes,'String'));
    temp.AzResUnits = get(handles.AzUnits,'String');
    temp.RangeResUnits = get(handles.RgUnits,'String');
    
    handles.recordparams = horzcat(handles.recordparams,temp);      
end

%if deployed, make sure image is updated now (do not want to do this
%normally, since it takes time)
if (isdeployed())
    drawnow expose;
end

% Update handles structure
guidata(handles.image, handles);