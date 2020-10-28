function varargout = ApertureTool(varargin)
%ApertureTool Interactive Image Aperture Viewer
%
% Interactive aperture tool.  Allows user to interactivly select the
% portion of the aperture that contributes to the image.  
%
% INPUTS:
%   filename   - optional : complex image filename
%   aoi        - optional : image AOI (filename must be specified)
%   segment    - optional : image number of multi-image file (filename must be specified)
%
% OUTPUTS:
%   GUI display, Movie and images as specified
%
% VERSION:
%   1.0
%     - Tim Cox 20100324
%     - Based on previous ApertureTool written by Tim Cox.
%       Initial concept based on tool written by Matt Banta, code lost to
%       antiquity
%   1.1
%     - Tim Cox 20110401
%     - Redesigned GUI.  Added resolution mode tab and functionaity.
%   1.2
%     - Wade Schwartzkopf 20121116
%     - Metadata-driven data normalization and general code cleanup
%   1.3
%     - Pat Cutler 20140124
%     - Added compatibility with polarized data
%     - Added image representation dropdown
%     - Added PHD axes labels
%   2.0 
%     - Tim Cox 20200909
%     - Got rid of old TabControl, uses new tab trick
%     - Redesigned tabs to move more common controls to first tab and
%       consolidated.  Moved some controls to menu items.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ApertureTool_OpeningFcn, ...
                   'gui_OutputFcn',  @ApertureTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ApertureTool is made visible.
function ApertureTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ApertureTool (see VARARGIN)

% Choose default command line output for ApertureTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.MaxLoadSize = 3000;

% Setup image and PHD axes
remap_list = getremaplist();
default_remap = find(strcmpi(remap_list,'densityremap'),1,'first');
set(handles.RemapCombo,'String',remap_list);
set(handles.RemapCombo,'Value',default_remap);
set(handles.ImageContrast,'Visible','off');
image_rep_list = {'amplitude' 'Pauli' 'AlphaEntropy' 'Stokes'};
set(handles.image_rep_combo,'String',image_rep_list);
set(handles.image_rep_combo,'Value',1);
set(handles.image_rep_combo,'Visible','off');
old_img_position = getpixelposition(handles.image);
handles.imghandle = imshow(0, 'Parent', handles.image);
handles.hSP1 = imscrollpanel(handles.figure1,handles.imghandle); 
setpixelposition(handles.hSP1,old_img_position);
handles.apiSP1 = iptgetapi(handles.hSP1);
set(handles.phd,'Visible','off');
set(handles.PHDPolCombo,'Visible','off');

%metaicon
set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');

%animation tab
set(handles.Slow,'Value',1);
set(handles.NumFrames,'String',7);
set(handles.ApFraction,'String',0.25);
set(handles.FrameRate,'String',5);
set(handles.CycleCheck,'Value',1);
set([handles.First,handles.Prev,...
    handles.Play,handles.Pause,...
    handles.Next,handles.Last],'enable','off');

%processing tab
set(handles.nx,'Enable','off');
set(handles.ny,'Enable','off');
set(handles.AlwaysApply,'Value',1);
set(handles.AlwaysApply,'Visible','off');


%set image/movie controls
set(handles.MovieFrameRate,'String',5);
set(handles.AnnotationCheck,'Value',1);
set(handles.StopRecording,'Enable','off');
set(handles.Recording,'Enable','off');
set(handles.SaveImage,'Enable','off');
set(handles.CurRecFrame,'Enable','off');
set(handles.TotRecFrame,'Enable','off');
set(handles.MovieFrameRate,'String',5);

%set up tab control
%handles.tgroup = uitabgroup('Parent', handles.figure1,'Position', [.023 .025 .497 .271]); 
handles.tgroup = uitabgroup('Parent', handles.figure1,'Position', [.023 .01 .51 .32]); 
handles.tab1 = uitab('Parent', handles.tgroup, 'Title', 'Animation');
handles.tab2 = uitab('Parent', handles.tgroup, 'Title', 'Processing');
handles.tab3 = uitab('Parent', handles.tgroup, 'Title', 'Save Image/Movie');
%handles.tab4 = uitab('Parent', handles.tgroup, 'Title', 'Resolution Mode');

%Place panels into each tab
set(handles.P1,'Parent',handles.tab1)
set(handles.P2,'Parent',handles.tab2)
%set(handles.P3,'Parent',handles.tab4)
set(handles.P4,'Parent',handles.tab3)

%Reposition each panel to same location as panel 1
set(handles.P2,'position',get(handles.P1,'position'));
%set(handles.P3,'position',get(handles.P1,'position'));
set(handles.P4,'position',get(handles.P1,'position'));

%add radio-button callbacks
set(handles.FilterPanel,'SelectionChangeFcn',@FilterChange);

%set status flags
handles.ResMode = 0;
handles.Record = 0;

% Update handles structure
guidata(hObject, handles);

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[0 0 0 0]);
p.addParamValue('segment',1);
p.parse(varargin{:});

%determine if Image name was passed as calling argument
if ~isempty(p.Results.filename)
    if ischar(p.Results.aoi) % Allows this to work deployed
        AOI = str2num(p.Results.aoi); % STR2NUM required since this will be a vector (must be passed in quotes for a deployed application)
    else
        AOI = p.Results.aoi;
    end
    if ischar(p.Results.segment)
        segment = str2double(p.Results.segment);
    else
        segment = p.Results.segment;
    end
    ApToolLoadImage(p.Results.filename,hObject,handles,AOI,segment);
end

% UIWAIT makes ApertureTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function FilterChange(hObject, eventdata)

handles = guidata(gcf);
ApToolnewpos(getPosition(handles.H), hObject);

% --- Outputs from this function are returned to the command line.
function varargout = ApertureTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PHDContrast.
function PHDContrast_Callback(hObject, eventdata, handles)
% hObject    handle to PHDContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imcontrast(handles.phd);

% --- Executes on selection change in PHDZoomCombo.
function PHDZoomCombo_Callback(hObject, eventdata, handles)
% hObject    handle to PHDZoomCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PHDZoomCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PHDZoomCombo
val = get(hObject,'Value');
% apiSP2 = iptgetapi(handles.hSP2);

%update image zoom
if (val == 1)
    set(handles.phd,'xlim',[1 handles.phdmagSize(2)],'ylim',[1 handles.phdmagSize(1)])
    %     set(handles.phdhandle,'cdata',handles.phdmag)
    %     mag = apiSP2.findFitMag();
    %     mag = mag-mag*0.2; %adjust for axes labels
    %     apiSP2.setMagnification(mag);
elseif (val == 2)
    xLims = handles.ZpLimsRn(1)+[0 handles.ZpWidthRn] ...
    +[-1 1]*handles.ZpWidthRn*.05; %add buffer
    xLims = [max(1,floor(xLims(1))) ...
        min(handles.phdmagSize(1),ceil(xLims(2)))]; %check image size
    yLims = handles.ZpLimsAz(1)+[0 handles.ZpWidthAz] ...
    +[-1 1]*handles.ZpWidthAz*.05; %add buffer
    yLims = [max(1,floor(yLims(1))) ...
        min(handles.phdmagSize(2),ceil(yLims(2)))]; %check image size
    set(handles.phd,'xlim',yLims,'ylim',xLims)
end
ApToolupdatePHD_axesLabels(handles)

% --- Executes during object creation, after setting all properties.
function PHDZoomCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PHDZoomCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RemapCombo.
function RemapCombo_Callback(hObject, eventdata, handles)
% hObject    handle to RemapCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

strings = get(handles.RemapCombo,'String');
val = get(handles.RemapCombo,'Value');

%make contrast button visible if remap is linearremap
if strcmpi(strings{val},'linearremap')
    set(handles.ImageContrast,'Visible','on');
else
    set(handles.ImageContrast,'Visible','off');
end

ApToolnewpos(getPosition(handles.H),hObject);

% --- Executes during object creation, after setting all properties.
function RemapCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RemapCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImageZoomCombo.
function ImageZoomCombo_Callback(hObject, eventdata, handles)
% hObject    handle to ImageZoomCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImageZoomCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImageZoomCombo
val = get(hObject,'Value');
apiSP1 = iptgetapi(handles.hSP1);

%update image zoom
if (val == 1)
    apiSP1.setMagnification(apiSP1.findFitMag());
else
    contents = cellstr(get(hObject,'String'));
    magPercent =  str2double(contents{val}(1:(end-1))); % 'end-1' removes percentage sign
    apiSP1.setMagnification(magPercent/100);
end


% --- Executes during object creation, after setting all properties.
function ImageZoomCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageZoomCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispref('matlab_sar_toolbox','last_used_directory')
    pathname = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathname)||~exist(pathname,'dir')
        pathname = pwd;
    end
else
    pathname = pwd;
end
[filenames,pathname]=uigetfile(sar_file_extensions('amplitude'),...
    'Open SAR Data File',pathname,'MultiSelect','on');

if ~isempty(pathname)
    setpref('matlab_sar_toolbox','last_used_directory',pathname); %store path
end

if(iscell(filenames)) % Multiple files requested
    for j=1:length(filenames)
        fullfilename{j}=[pathname filenames{j}];
    end
elseif(filenames)
    fullfilename{1}=[pathname filenames];
end

ApToolLoadImage(fullfilename,hObject,handles,[],1);

% --- Executes on button press in Measure.
function Measure_Callback(hObject, eventdata, handles)
% hObject    handle to Measure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject,'String'),'Measure')
    handles.MeasureLine = imline(handles.image);
    setColor(handles.MeasureLine,'Cyan');
    kids = get(handles.MeasureLine,'Children');
    set(kids,'LineWidth',3);
    set(kids,'MarkerSize',1);
    api_handle = iptgetapi(handles.MeasureLine);
    api_handle.addApToolnewpositionCallback(@newlinepos);
    set(handles.Measure,'String','Remove');
    %compute line distance and heading in the ground plane
    pos = getPosition(handles.MeasureLine);
    StartPos = pos(1,:);
    StopPos = pos(2,:);
    StartPos(1) = StartPos(1)+handles.aoi(1);
    StartPos(2) = StartPos(2)+handles.aoi(2);
    StopPos(1) = StopPos(1)+handles.aoi(1);
    StopPos(2) = StopPos(2)+handles.aoi(2);
    lla1 = point_slant_to_ground(fliplr(StartPos)',handles.meta);
    lla2 = point_slant_to_ground(fliplr(StopPos)',handles.meta);
    d1km=lldistkm(lla1(1:2),lla2(1:2));
    Distance = d1km*1000; %meters
    %get units (meters for metric and feet for English
    if get(handles.MetricCheck,'Value')
        set(handles.MeasureDist,'String',sprintf('%.1f',Distance));
        set(handles.MeasureUnits,'String','meters');        
    else
        set(handles.MeasureDist,'String',sprintf('%.1f',Distance/0.3048));
        set(handles.MeasureUnits,'String','feet'); 
    end
    set(handles.MeasureUnits,'Visible','on');        
else
    delete(handles.MeasureLine);
    handles.MeasureLine = [];
    set(handles.Measure,'String','Measure');
    set(handles.MeasureDist,'String','');
    set(handles.MeasureUnits,'Visible','on'); 
end

guidata(hObject, handles);

function MeasureDist_Callback(hObject, eventdata, handles)
% hObject    handle to MeasureDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MeasureDist as text
%        str2double(get(hObject,'String')) returns contents of MeasureDist as a double


% --- Executes during object creation, after setting all properties.
function MeasureDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MeasureDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FullAperture.
function FullAperture_Callback(hObject, eventdata, handles)
% hObject    handle to FullAperture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos = [handles.ZpLimsAz(1), handles.ZpLimsRn(1),...
    handles.ZpWidthAz handles.ZpWidthRn];
setPosition(handles.H,pos);
ApToolnewpos(pos, hObject);

function AzRes_Callback(hObject, eventdata, handles)
% hObject    handle to AzRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AzRes as text
%        str2double(get(hObject,'String')) returns contents of AzRes as a double


% --- Executes during object creation, after setting all properties.
function AzRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AzRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RgRes_Callback(hObject, eventdata, handles)
% hObject    handle to RgRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RgRes as text
%        str2double(get(hObject,'String')) returns contents of RgRes as a double


% --- Executes during object creation, after setting all properties.
function RgRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RgRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MetricCheck.
function MetricCheck_Callback(hObject, eventdata, handles)
% hObject    handle to MetricCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MetricCheck

ApToolnewpos(getPosition(handles.H),hObject);

function CurFrame_Callback(hObject, eventdata, handles)
% hObject    handle to CurFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurFrame as text
%        str2double(get(hObject,'String')) returns contents of CurFrame as a double


% --- Executes during object creation, after setting all properties.
function CurFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TotFrames_Callback(hObject, eventdata, handles)
% hObject    handle to TotFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TotFrames as text
%        str2double(get(hObject,'String')) returns contents of TotFrames as a double


% --- Executes during object creation, after setting all properties.
function TotFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in First.
function First_Callback(hObject, eventdata, handles)
% hObject    handle to First (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set([handles.Next handles.Last],'enable','on');
set([handles.First handles.Prev],'enable','off');
UpdateAnimationFrame(1, handles);

% --- Executes on button press in Prev.
function Prev_Callback(hObject, eventdata, handles)
% hObject    handle to Prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%previous Frame
CurFrame = str2double(get(handles.CurFrame,'String')) - 1;

set([handles.Next handles.Last],'enable','on');

if (CurFrame == 1)    
    set([handles.First handles.Prev],'enable','off');
else
    set([handles.First handles.Prev],'enable','on');
end

UpdateAnimationFrame(CurFrame, handles);

% --- Executes on button press in Play.
function Play_Callback(hObject, eventdata, handles)
% hObject    handle to Play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get settings
NumFrames = str2double(get(handles.NumFrames,'String'));
FrameRate = str2double(get(handles.MovieFrameRate,'String'));
CycleCheck = get(handles.CycleCheck,'Value'); % Cycling continuously?

%disable all controls except pause
set(handles.First,'enable','off');
set(handles.Prev,'enable','off');
set(handles.Play,'enable','off');
set(handles.Pause,'enable','on');
set(handles.Next,'enable','off');
set(handles.Last,'enable','off');

% Set flag we will use for recording
handles.InLoop = 1;
guidata(hObject, handles);

set(handles.TotFrames,'String',NumFrames);
CurFrame = 1;

while (handles.InLoop)
    
    %see if pause was pressed
    if (strcmp(get(handles.Play,'Enable'),'on'))
        handles.InLoop = 0;
        break;
    end
    
    tic;
    UpdateAnimationFrame(CurFrame, handles);
    %wait until time for next frame
    elapsedtime = toc;
    pausetime = 1/FrameRate - elapsedtime;                
    if (pausetime > 0)
        pause(pausetime);
    else        
        pause(.05); %pause a little so the animation can be stopped
    end
    
    if (CurFrame == NumFrames)
        if (CycleCheck)
            CurFrame = 1;
        else
            %stop movie
            set(handles.Play,'enable','on');
            set(handles.Pause,'enable','off');
            break;
        end
    else
        CurFrame = CurFrame + 1;
    end
end

% --- Executes on button press in Pause.
function Pause_Callback(hObject, eventdata, handles)
% hObject    handle to Pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Play,'Enable','on');
set(handles.NumFrames,'Enable','on');
set(handles.ApFraction,'Enable','on');
set(handles.MovieFrameRate,'Enable','on');

CurFrame = str2double(get(handles.CurFrame,'String'));
TotFrames = str2double(get(handles.TotFrames,'String'));

if (CurFrame > 1)
    set([handles.Prev handles.First],'Enable','on');
end

if (CurFrame < TotFrames)
    set([handles.Next handles.Last],'Enable','on');
end

set(handles.Pause,'Enable','off');

% --- Executes on button press in Next.
function Next_Callback(hObject, eventdata, handles)
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%advance Frame
CurFrame = str2double(get(handles.CurFrame,'String')) + 1;
TotFrames = str2double(get(handles.NumFrames,'String'));

if (CurFrame > 1)
    set([handles.First handles.Prev],'enable','on');
else
    set([handles.First handles.Prev],'enable','off');
end

if (CurFrame == TotFrames)    
    set([handles.Next handles.Last],'enable','off');
else
    set([handles.Next handles.Last],'enable','on');
end

UpdateAnimationFrame(CurFrame, handles);

% --- Executes on button press in Last.
function Last_Callback(hObject, eventdata, handles)
% hObject    handle to Last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%advance Frame
set([handles.First handles.Prev],'enable','on');
set([handles.Next handles.Last],'enable','off');

CurFrame = str2double(get(handles.NumFrames,'String'));
UpdateAnimationFrame(CurFrame, handles);

function UpdateAnimationFrame(CurFrame, handles)

if isnan(CurFrame)
    CurFrame = 1;
end
set(handles.CurFrame,'String',CurFrame);

NumFrames = str2double(get(handles.NumFrames,'String'));
ApertureFraction = str2double(get(handles.ApFraction,'String'));
if (get(handles.Slow,'Value'))
    ApWidth = round(handles.ZpWidthAz*ApertureFraction);
    StartPos = max(1,handles.ZpLimsAz(1) + ...
        round((CurFrame-1)*(handles.ZpWidthAz-ApWidth)/(NumFrames-1)));

    %update imrect (this will update image)
    setPosition(handles.H,[StartPos handles.ZpLimsRn(1) ApWidth handles.ZpWidthRn]);
else
    ApWidth = round(handles.ZpWidthRn*ApertureFraction);
    StartPos = max(1,handles.ZpLimsRn(1) + ...
        round((CurFrame-1)*(handles.ZpWidthRn-ApWidth)/(NumFrames-1)));
    [ny, nx] = size(handles.phasehistory);
    StopPos = min(ny,StartPos + ApWidth);
    StartPos = ny - StopPos + 1; % Flip due to imshow convention

    %update imrect        
    setPosition(handles.H,[handles.ZpLimsAz(1) StartPos handles.ZpWidthAz ApWidth]);        
end


function ApTime_Callback(hObject, eventdata, handles)
% hObject    handle to ApTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApTime as text
%        str2double(get(hObject,'String')) returns contents of ApTime as a double


% --- Executes during object creation, after setting all properties.
function ApTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ApBandwidth_Callback(hObject, eventdata, handles)
% hObject    handle to ApBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApBandwidth as text
%        str2double(get(hObject,'String')) returns contents of ApBandwidth as a double


% --- Executes during object creation, after setting all properties.
function ApBandwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DeskewSlow.
function DeskewSlow_Callback(hObject, eventdata, handles)
% hObject    handle to DeskewSlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') % In general, we can't deskew in range and cross-range simultaneously
    if ~isscalar(handles.meta.Grid.Col.DeltaKCOAPoly) || ...
            ~isscalar(handles.meta.Grid.Row.DeltaKCOAPoly)
        set(handles.DeskewFast,'Value',false);
    end
    handles.manual_offset = [0 0];
elseif all(handles.meta.Grid.Row.DeltaKCOAPoly(:)==0)&&handles.isvalid_row_deltakcoapoly
    set(handles.DeskewFast,'Value',true);
end
ApToolprocess_phd(hObject,handles);


% --- Executes on button press in DeskewFast.
function DeskewFast_Callback(hObject, eventdata, handles)
% hObject    handle to DeskewFast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') % In general, we can't deskew in range and cross-range simultaneously
    if ~isscalar(handles.meta.Grid.Col.DeltaKCOAPoly) || ...
            ~isscalar(handles.meta.Grid.Row.DeltaKCOAPoly)
        set(handles.DeskewSlow,'Value',false);
    end
    handles.manual_offset = [0 0];
elseif all(handles.meta.Grid.Col.DeltaKCOAPoly(:)==0) && handles.isvalid_col_deltakcoapoly
    set(handles.DeskewSlow,'Value',true);
end
ApToolprocess_phd(hObject,handles);


% --- Executes on button press in UniformWeighting.
function UniformWeighting_Callback(hObject, eventdata, handles)
% hObject    handle to UniformWeighting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UniformWeighting
ApToolprocess_phd(hObject,handles);

% --- Executes on button press in ManualAdjust.
function ManualAdjust_Callback(hObject, eventdata, handles)
% hObject    handle to ManualAdjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fhand = figure('Name','Select Nonzeropad Area','Colormap',gray,...
    'MenuBar','none','Toolbar','none','NumberTitle','off');
ahand = axes('Parent',fhand,'Units','normalized','Position',[0 .2 1 .8]);
uicontrol('Style','pushbutton','String','Recenter','FontWeight','bold',...
    'Units','normalized','Position',[.2 .05 .2 .1],'parent',fhand,...
    'Callback',@manual_recenter);
uicontrol('Style','pushbutton','String','OK','FontWeight','bold',...
    'Units','normalized','Position',[.6 .05 .2 .1],'parent',fhand,...
    'Callback',@manual_ok);
if isfield(handles,'polar_PHD_idx')
    local.phasehistory = abs(fftshift(handles.fft_sp(handles.complex_data(:,:,handles.polar_PHD_idx).')));
else
    local.phasehistory = abs(fftshift(handles.fft_sp(handles.complex_data.')));
end
% PHD display is flipped so that time increase to the right
if isfield(handles.meta,'SCPCOA') && isfield(handles.meta.SCPCOA,'SideOfTrack') && ...
        strcmp(handles.meta.SCPCOA.SideOfTrack,'L')
    local.phasehistory = fliplr(local.phasehistory);
end

% We could guess frequency support center/zeropad from data here to
% intialize imrect position.
phdmean = mean(local.phasehistory(:));
phdstd = std(local.phasehistory(:));
local.ihand = image(local.phasehistory,'Parent',ahand,'CDataMapping','scaled');
set(ahand,'DataAspectRatio',[1 1 1],'CLim',[0 phdmean+3*phdstd],...
    'XTick',[],'YTick',[]);
local.manual_imrect = imrect(ahand, [handles.ZpLimsAz(1), handles.ZpLimsRn(1),...
    handles.ZpWidthAz handles.ZpWidthRn]);
local.manual_offset = [0 0]; % Remember how offcenter we are
local.parent = hObject;
guidata(fhand, local);


% Recenter image of PHD to imrect
function manual_recenter(hObject, eventdata)
handles = guidata(hObject);
[ny, nx] = size(handles.phasehistory);
center = round([nx,ny]/2);
pos = getPosition(handles.manual_imrect);
if all(pos([1 2])<=1) && all((pos([1 2])+pos([3 4]))>[nx ny])
    pos = [1 1 nx ny];
end
handles.manual_offset = mod(center - round(pos([1 2])+(pos([3 4])/2)) + ...
    + handles.manual_offset,[nx ny]);
guidata(hObject, handles); % Update handles.manual_offset
set(handles.ihand,'CData',circshift(handles.phasehistory,handles.manual_offset([2 1])));
pos([1 2]) = center - round(pos([3 4])/2);
setPosition(handles.manual_imrect,pos);


function manual_ok(hObject, eventdata)
manual_recenter(hObject, eventdata);
handles = guidata(hObject);
parent = guidata(handles.parent);
[ny, nx] = size(handles.phasehistory);
pos = getPosition(handles.manual_imrect);
parent.AzPad = max(1,nx/pos(3));
parent.RnPad = max(1,ny/pos(4));
parent.ZpLimsAz = max(1,pos(1));
parent.ZpLimsAz(2) = min(nx,pos(1)+pos(3)-1);
parent.ZpWidthAz = diff(parent.ZpLimsAz)+1;
parent.ZpLimsRn = max(1,pos(2));
parent.ZpLimsRn(2) = min(ny,pos(2)+pos(4)-1);
parent.ZpWidthRn = diff(parent.ZpLimsRn)+1;

if verLessThan('matlab', '7.6')
    disp('Warning. Matlab version does not support imrect constraint function');
else
    setPositionConstraintFcn(parent.H,...
        makeConstrainToRectFcn('imrect',parent.ZpLimsAz,parent.ZpLimsRn));
end

parent.manual_offset = handles.manual_offset;
% PHD display is flipped so that time increase to the right, but offset is
% applied before flipping.
if isfield(parent.meta,'SCPCOA') && isfield(parent.meta.SCPCOA,'SideOfTrack') && ...
        strcmp(parent.meta.SCPCOA.SideOfTrack,'L')
    parent.manual_offset(1) = -parent.manual_offset(1);
end

% Recompute weighting here.  If frequency support/zeropad factor were wrong
% (which they must have been, if the data has to be manually adjusted), and
% the weighting was estimated from the pixel data (as opposed to the
% metadata), the weighting function was almost certainly computed
% incorrectly.
if parent.weighting_computed_from_data
    phasehistory = parent.fft_sp(conj(parent.complex_data));
    phasehistory = circshift(phasehistory,handles.manual_offset);
    data = parent.fft_im(phasehistory);
    parent.weight_fun_az = estimate_weighting_mem(data, 1, parent.AzPad);
    parent.weight_fun_rng = estimate_weighting_mem(data, 2, parent.RnPad);
end

set(parent.DeskewSlow,'Value',false);
set(parent.DeskewFast,'Value',false);
ApToolprocess_phd(handles.parent,parent);
close(get(hObject,'Parent'));

% --- Executes on button press in CycleCheck.
function CycleCheck_Callback(hObject, eventdata, handles)
% hObject    handle to CycleCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CycleCheck



function NumFrames_Callback(hObject, eventdata, handles)
% hObject    handle to NumFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumFrames as text
%        str2double(get(hObject,'String')) returns contents of NumFrames as a double


% --- Executes during object creation, after setting all properties.
function NumFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ApFraction_Callback(hObject, eventdata, handles)
% hObject    handle to ApFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApFraction as text
%        str2double(get(hObject,'String')) returns contents of ApFraction as a double


% --- Executes during object creation, after setting all properties.
function ApFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MovieFrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to MovieFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MovieFrameRate as text
%        str2double(get(hObject,'String')) returns contents of MovieFrameRate as a double


% --- Executes during object creation, after setting all properties.
function MovieFrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MovieFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ImageContrast.
function ImageContrast_Callback(hObject, eventdata, handles)
% hObject    handle to ImageContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imcontrast(handles.image);

% --- Executes on selection change in image_rep_combo.
function image_rep_combo_Callback(hObject, eventdata, handles)
% hObject    handle to image_rep_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_rep_combo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_rep_combo
ApToolnewpos(getPosition(handles.H),hObject);

% --- Executes during object creation, after setting all properties.
function image_rep_combo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_rep_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PHDPolCombo.
function PHDPolCombo_Callback(hObject, eventdata, handles)
% hObject    handle to PHDPolCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PHDPolCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PHDPolCombo

handles.polar_PHD_idx = get(handles.PHDPolCombo,'Value');
ApToolprocess_phd(hObject,handles); % Apply PHD options

% --- Executes during object creation, after setting all properties.
function PHDPolCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PHDPolCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in InverseSelection.
function InverseSelection_Callback(hObject, eventdata, handles)
% hObject    handle to InverseSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InverseSelection
if get(hObject,'Value')
    setColor(handles.H,'red');
else
    setColor(handles.H,'blue');
end
ApToolnewpos(getPosition(handles.H),hObject);

% --- Executes on button press in InversePolar.
function InversePolar_Callback(hObject, eventdata, handles)
% hObject    handle to InversePolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InversePolar
ApToolprocess_phd(hObject,handles); % Apply PHD options
ApToolnewpos(getPosition(handles.H),hObject);

% --- Executes on button press in AlwaysApply.
function AlwaysApply_Callback(hObject, eventdata, handles)
% hObject    handle to AlwaysApply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AlwaysApply
ApToolnewpos(getPosition(handles.H),hObject);


function nx_Callback(hObject, eventdata, handles)
% hObject    handle to nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nx as text
%        str2double(get(hObject,'String')) returns contents of nx as a double


% --- Executes during object creation, after setting all properties.
function nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ny_Callback(hObject, eventdata, handles)
% hObject    handle to ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ny as text
%        str2double(get(hObject,'String')) returns contents of ny as a double


% --- Executes during object creation, after setting all properties.
function ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AnnotationCheck.
function AnnotationCheck_Callback(hObject, eventdata, handles)
% hObject    handle to AnnotationCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AnnotationCheck



function CurRecFrame_Callback(hObject, eventdata, handles)
% hObject    handle to CurRecFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurRecFrame as text
%        str2double(get(hObject,'String')) returns contents of CurRecFrame as a double


% --- Executes during object creation, after setting all properties.
function CurRecFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurRecFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TotRecFrame_Callback(hObject, eventdata, handles)
% hObject    handle to TotRecFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TotRecFrame as text
%        str2double(get(hObject,'String')) returns contents of TotRecFrame as a double


% --- Executes during object creation, after setting all properties.
function TotRecFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotRecFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Recording.
function Recording_Callback(hObject, eventdata, handles)
% hObject    handle to Recording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Record = 1;
set(handles.Recording,'enable','off');
set(handles.StopRecording,'enable','on');
handles.recordparams = [];

%turn of cycle continuously check as this is not desired for a movie (the
%GIF will cycle on its' own)
set(handles.CycleCheck,'Value',0);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in StopRecording.
function StopRecording_Callback(hObject, eventdata, handles)
% hObject    handle to StopRecording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get selection and launch appropriate file dialog
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

filterspec={'*.gif','Animated GIF (*.gif)'};
[fname, path] = uiputfile( filterspec,'Save Movie File',pathstr);
if isequal(fname,0) % Cancel was pressed
    return;
end

filename = strcat(path,fname);
[pathstr, name, ext] = fileparts(filename);

setpref('matlab_sar_toolbox','last_used_directory',path); %store path

set(handles.Recording,'Enable','off');
set(handles.StopRecording,'Enable','off');

numframes = length(handles.recordparams);
set(handles.TotRecFrame,'String',numframes);

%get settings
FullRes = get(handles.FullRes,'Value');
Annotate = get(handles.AnnotationCheck,'Value');
framerate = str2double(get(handles.MovieFrameRate,'String'));

if FullRes
    [chip, handles] = getFullResImage(handles);
    [ny,nx,~] = size(chip);

    %if chip is small then we will make bigger if annotation is selected,
    %otherwise it will not be legible
    MinSize = 1000;
    if Annotate
        if nx<MinSize && ny<MinSize
            if ny < nx
                new_ny = MinSize;
                new_nx = round((new_ny/ny)*nx);
            else
                new_nx = MinSize;
                new_ny = round((new_nx/nx)*ny);
            end
            chip = imresize(chip,[new_ny new_nx]);
        end
    end
else % Screen resolution doesn't do annotation
    F = getframe( handles.image );
    chip = F.cdata(:,:,:);
end

val = get(handles.image_rep_combo,'Value');
strings = get(handles.image_rep_combo,'String');
polstring = strings{val};
animated(:,:,:,numframes) = chip;

TextCount = 0;

for i=1:numframes
    if FullRes
        %Full Resolution, get full chip for
        setPosition(handles.H,handles.recordparams(i).Pos);
        chip = abs(getFullResImage(handles));
        
        if Annotate
            %Compute Height of Text based on ny, we'll use 1/8th
            %Mode
            %Frame n of n
            %AzRes: # units, RgRes: # units
            BlockHeight = round(size(chip,1)/8); 
            FontSize = round(BlockHeight/3);
            ModeString = ['Mode: ' handles.recordparams(i).Mode];
            ModeStringIm = RasterizeText(ModeString,'Times New Roman',FontSize);
            FrameString = ['Frame ' num2str(i) ' of ' num2str(numframes)];
            FrameStringIm = RasterizeText(FrameString,'Times New Roman',FontSize);
            ResString = sprintf('AzRes: %3.1f %s, RgRes: %3.1f %s', handles.recordparams(i).AzRes, ...
                        handles.recordparams(i).AzResUnits, handles.recordparams(i).RangeRes, ...
                        handles.recordparams(i).RangeResUnits);
            ResStringIm = RasterizeText(ResString,'Times New Roman',FontSize);
            ny = size(ModeStringIm,1)+size(FrameStringIm,1)+size(ResStringIm,1);
            nx = max([size(ModeStringIm,2) size(FrameStringIm,2) size(ResStringIm,2)]);
            TextIm = zeros(ny,nx);
            TextIm(1:size(ModeStringIm,1),1:size(ModeStringIm,2)) = ModeStringIm;
            StartRow = size(ModeStringIm,1)+1;
            StopRow = StartRow + size(FrameStringIm,1)-1;
            TextIm(StartRow:StopRow,1:size(FrameStringIm,2)) = FrameStringIm;
            StartRow = StopRow+1;
            StopRow = StartRow + size(ResStringIm,1)-1;
            TextIm(StartRow:StopRow,1:size(ResStringIm,2)) = ResStringIm;
            TextIm(TextIm>0)=255;
            %Add annotation to lower left
            StartRow = size(chip,1)-ny+1;
            if size(chip,3) == 1
                chip(StartRow:end,1:nx) = TextIm;
            else
                chip(StartRow:end,1:nx,1) = TextIm;
                chip(StartRow:end,1:nx,2) = TextIm;
                chip(StartRow:end,1:nx,3) = TextIm;
            end
        end
    else
        %Screen Resolution, just update imrect pos and let callback update
        setPosition(handles.H,handles.recordparams(i).Pos);
        %set CLim
        set(handles.image,'CLim',handles.recordparams(i).CLim);
        %set scroll panel to selected view
        handles.apiSP1.setMagnificationAndCenter(handles.recordparams(i).Mag,...
        handles.recordparams(i).Center(1),handles.recordparams(i).Center(2));
        axes(handles.image);
        F = getframe( handles.image );
        chip = F.cdata(:,:,:);

        if Annotate
            %Compute Height of Text based on ny, we'll use 1/8th
            %Mode
            %Frame n of n
            %AzRes: # units, RgRes: # units
            BlockHeight = round(size(chip,1)/8); 
            FontSize = round(BlockHeight/3);
            ModeString = ['Mode: ' handles.recordparams(i).Mode];
            ModeStringIm = RasterizeText(ModeString,'Times New Roman',FontSize);
            FrameString = ['Frame ' num2str(i) ' of ' num2str(numframes)];
            FrameStringIm = RasterizeText(FrameString,'Times New Roman',FontSize);
            ResString = sprintf('AzRes: %3.1f %s, RgRes: %3.1f %s', handles.recordparams(i).AzRes, ...
                        handles.recordparams(i).AzResUnits, handles.recordparams(i).RangeRes, ...
                        handles.recordparams(i).RangeResUnits);
            ResStringIm = RasterizeText(ResString,'Times New Roman',FontSize);
            ny = size(ModeStringIm,1)+size(FrameStringIm,1)+size(ResStringIm,1);
            nx = max([size(ModeStringIm,2) size(FrameStringIm,2) size(ResStringIm,2)]);
            TextIm = zeros(ny,nx);
            TextIm(1:size(ModeStringIm,1),1:size(ModeStringIm,2)) = ModeStringIm;
            StartRow = size(ModeStringIm,1)+1;
            StopRow = StartRow + size(FrameStringIm,1)-1;
            TextIm(StartRow:StopRow,1:size(FrameStringIm,2)) = FrameStringIm;
            StartRow = StopRow+1;
            StopRow = StartRow + size(ResStringIm,1)-1;
            TextIm(StartRow:StopRow,1:size(ResStringIm,2)) = ResStringIm;
            TextIm(TextIm>0)=255;
            %Add annotation to lower left
            StartRow = size(chip,1)-ny+1;
            chip(StartRow:end,1:nx,1) = TextIm;
            chip(StartRow:end,1:nx,2) = TextIm;
            chip(StartRow:end,1:nx,3) = TextIm;
        end
    end
       
    if ~strcmpi(polstring,'amplitude')
        %convert to indexed image
        if i == 1
            [animated,cmap] = rgb2ind(chip,256,'nodither');
        else
            animated(:,:,1,i) = rgb2ind(chip,cmap,'nodither');
        end
    else
        if i == 1
            animated(:,:,:,i) = double(chip);
        else
            animated(:,:,:,i) = double(chip);
        end
    end
    
    %update status
    set(handles.CurRecFrame,'String',i);    
end

DelayTime = 1/framerate;
if strcmpi(polstring,'amplitude')
    imwrite(animated, filename, 'DelayTime', DelayTime, 'LoopCount', inf);
else
    imwrite(animated, cmap, filename, 'DelayTime', DelayTime, 'LoopCount', inf);
end

set(handles.TotRecFrame,'String','');
set(handles.CurRecFrame,'String','');

set(handles.Recording,'Enable','on');
set(handles.StopRecording,'Enable','off');

handles.Record = 0;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SaveImage.
function SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

if get(handles.JPG,'Value')
    [fname, path] = uiputfile( {'*.jpg','jpg Files (*.jpg)'},'Save JPG File',pathstr);
    filename = strcat(path,fname);
    if get(handles.FullRes,'Value')
        %full res image of chip
        chip = abs(getFullResImage(handles));
        imwrite(chip,filename,'jpg');
    else
        %just write what is on the screen
        frame = getframe(handles.image);
        imwrite(frame.cdata,filename,'jpg');
    end   
else
    [fname, path] = uiputfile( {'*.ntf','SICD Files (*.ntf)'},'Save SICD File',pathstr);
    filename = strcat(path,fname);    

    %get complex image chip
    chip = getFullResImage(handles,'complex');
    chip = handles.fft_sp(fftshift(handles.fft_im((chip))));
    meta = handles.meta;
    meta.ImageFormation.Processing.Type = 'Subaperture';

    %update res, frequency support etc.
    %this should work for spotlight collects

    %get new BW and center position
    pos = getPosition(handles.H);
    SlowPercent = pos(3)/handles.ZpWidthAz;
    SlowCenter = (pos(1)+pos(3)/2)/size(chip,2);
    FastPercent = pos(4)/handles.ZpWidthRn;
    FastCenter = (pos(2)+pos(4)/2)/size(chip,1);

    %update SICD fields
    meta.Grid.Col.ImpRespBW = meta.Grid.Col.ImpRespBW*SlowPercent;
    meta.Grid.Col.DeltaKCOAPoly = double(meta.Grid.Col.Sgn)*(SlowCenter-.5)/meta.Grid.Col.SS;
    meta.Grid.Row.ImpRespBW = meta.Grid.Row.ImpRespBW*FastPercent;
    meta.Grid.Row.DeltaKCOAPoly = double(meta.Grid.Row.Sgn)*(FastCenter-.5)/meta.Grid.Row.SS;
    %set remaining fields to zero and re-derive
    meta.Grid.Row = rmfield(meta.Grid.Row,'ImpRespWid');
    meta.Grid.Row = rmfield(meta.Grid.Row,'DeltaK1');
    meta.Grid.Row = rmfield(meta.Grid.Row,'DeltaK2');
    meta.Grid.Col = rmfield(meta.Grid.Col,'ImpRespWid');
    meta.Grid.Col = rmfield(meta.Grid.Col,'DeltaK1');
    meta.Grid.Col = rmfield(meta.Grid.Col,'DeltaK2');
    meta = derived_sicd_fields(meta);

    %if weighting was applied, then set the field...TODO: add sampled
    %weighting function so it can be undone
    %we will assume any weighting from the input image has been removed
    if get(handles.Gaussian,'Value')
        WgtType = 'GAUSSIAN';
    elseif get(handles.x4,'Value')
        WgtType = '1/X4';
    elseif get(handles.Hamming,'Value')
        WgtType = 'HAMMING';
    elseif get(handles.CosOnPed,'Value')
        WgtType = 'COSINEONPED';
    else
        WgtType = 'UNIFORM';
    end
    meta.Grid.Col.WgtType = WgtType;
    meta.Grid.Row.WgtType = WgtType;

    meta.ImageData.NumRows = size(chip,1);
    meta.ImageData.NumCols = size(chip,2);
    meta = add_sicd_corners(meta); 

    writer = SICDWriter(filename,meta);
    writer.write_chip(chip.',[1 1]);
    delete(writer);    
end

function [chip,handles] = getFullResImage(handles,format)

if ~exist('format','var'); format = 'detected'; end

%return full res complex image chip for current position
pos = getPosition(handles.H);
xmin = max(handles.ZpLimsAz(1),floor(pos(1)));
xmax = min(handles.ZpLimsAz(2),floor(pos(1)+pos(3)-1));
ymin = max(handles.ZpLimsRn(1),floor(pos(2)));
ymax = min(handles.ZpLimsRn(2),floor(pos(2)+pos(4)-1));

%need to modify xposition if flight is left
[ny,nx,~] = size(handles.phasehistory);
if isfield(handles.meta,'SCPCOA') && isfield(handles.meta.SCPCOA,'SideOfTrack') && ...
   strcmp(handles.meta.SCPCOA.SideOfTrack,'L')
    xminold = xmin;
    xmin = round(nx-xmax);
    if (xmin<1);xmin=1;end
    xmax = round(nx-xminold);    
end

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
    FilterRange = 1;
    FilterAz = 1;
end

%now reform image with selected aperture
if get(handles.InverseSelection,'Value')        
    Ap = single(handles.phasehistory);
    Ap(ymin:ymax,xmin:xmax) = 0;
else
    Ap = zeros(ny,nx,'single');
    phdchip = ApToolFilterPHDChip(phdchip,WindowFilter,FilterRange,FilterAz); 
    Ap(ymin:ymax,xmin:xmax,1:size(phdchip,3)) = phdchip;    
end

if ~strcmpi(format,'complex')
    [chip,handles] = ApToolmakeDisplayable(handles,Ap);
else
    chip = handles.fft_im(Ap);
end



function FrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to FrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameRate as text
%        str2double(get(hObject,'String')) returns contents of FrameRate as a double


% --- Executes during object creation, after setting all properties.
function FrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
