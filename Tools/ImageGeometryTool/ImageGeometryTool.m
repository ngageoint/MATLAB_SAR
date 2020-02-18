function varargout = ImageGeometryTool(varargin)
% IMAGEGEOMETRYTOOL MATLAB code for ImageGeometryTool.fig
%      IMAGEGEOMETRYTOOL, by itself, creates a new IMAGEGEOMETRYTOOL or raises the existing
%      singleton*.
%
%      H = IMAGEGEOMETRYTOOL returns the handle to a new IMAGEGEOMETRYTOOL or the handle to
%      the existing singleton*.
%
%      IMAGEGEOMETRYTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEGEOMETRYTOOL.M with the given input arguments.
%
%      IMAGEGEOMETRYTOOL('Property','Value',...) creates a new IMAGEGEOMETRYTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageGeometryTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageGeometryTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageGeometryTool

% Last Modified by GUIDE v2.5 07-Feb-2020 09:58:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageGeometryTool_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageGeometryTool_OutputFcn, ...
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


% --- Executes just before ImageGeometryTool is made visible.
function ImageGeometryTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageGeometryTool (see VARARGIN)

% Choose default command line output for ImageGeometryTool
handles.output = hObject;

set(handles.image,'XTickLabel',[]);
set(handles.image,'YTickLabel',[]);

set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');
set(handles.metaicon,'Color',[51/255 102/255 153/255]);

handles.ImageLoaded = 0;
handles.ChangingView = 0;
handles.buttondown = 0;

%turn everything on to start
set(handles.VelocityCheck,'Value',1);
set(handles.RangeCheck,'Value',1);
set(handles.RangeGroundCheck,'Value',1);
set(handles.SlantPlaneCheck,'Value',1);
set(handles.GroundPlaneCheck,'Value',1);
set(handles.GroundBorderCheck,'Value',1);
set(handles.GroundSlantCheck,'Value',1);
set(handles.SSPCheck,'Value',1);
set(handles.GroundNormalCheck,'Value',1);
set(handles.GroundLayoverCheck,'Value',1);
set(handles.GroundMultipathCheck,'Value',1);
set(handles.GroundShadowCheck,'Value',1);
set(handles.GroundNorthCheck,'Value',1);
set(handles.SlantNormalCheck,'Value',1);
set(handles.SlantLayoverCheck,'Value',1);
set(handles.SlantMultipathCheck,'Value',1);
set(handles.SlantShadowCheck,'Value',1);
set(handles.SlantNorthCheck,'Value',1);

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[0 0 0 0]);
p.addParamValue('segment',1);
p.parse(varargin{:});

if length(p.Results.filename) > 0
    if isfield(p.Results,'aoi')
        handles.AOI = p.Results.aoi;
    else
        handles.AOI = [0 0 0 0];
    end
    if isfield(p.Results,'segment')
        handles.segment =  p.Results.segment;
    else
        handles.segment = 1;
    end
    handles.filename = p.Results.filename;
    guidata(hObject, handles);
    LoadImage(handles);       
else
    % Update handles structure
    guidata(hObject, handles);
end

% UIWAIT makes ImageGeometryTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function LoadImage(handles)

reader_obj = open_reader(handles.filename);
if iscell(reader_obj)
    reader_obj = reader_obj{1};
end
handles.meta = reader_obj.get_meta();

try
    MetaIcon_Complex(handles.meta,'handle',handles.metaicon);
catch
end

UpdateView(handles);

function ViewSelect(hObject, eventdata)

handles = guidata(gcf);

%disable sub-aperture controls for Perspective View
if get(handles.PerspectiveView,'Value')
    set(handles.SubApPanel,'Visible','off');    
else
    set(handles.SubApPanel,'Visible','on');
end

UpdateView(handles);



% --- Outputs from this function are returned to the command line.
function varargout = ImageGeometryTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in VelocityCheck.
function VelocityCheck_Callback(hObject, eventdata, handles)
% hObject    handle to VelocityCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VelocityCheck
UpdateView(handles);

% --- Executes on button press in SlantPlaneCheck.
function SlantPlaneCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlantPlaneCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlantPlaneCheck
UpdateView(handles);

% --- Executes on button press in RangeCheck.
function RangeCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RangeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RangeCheck
UpdateView(handles);

% --- Executes on button press in RangeGroundCheck.
function RangeGroundCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RangeGroundCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RangeGroundCheck
UpdateView(handles);

% --- Executes on button press in GroundPlaneCheck.
function GroundPlaneCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundPlaneCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundPlaneCheck
UpdateView(handles);

% --- Executes on button press in GroundSlantCheck.
function GroundSlantCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundSlantCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundSlantCheck
UpdateView(handles);

% --- Executes on button press in GroundNormalCheck.
function GroundNormalCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundNormalCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundNormalCheck
UpdateView(handles);

% --- Executes on button press in GroundLayoverCheck.
function GroundLayoverCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundLayoverCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundLayoverCheck
UpdateView(handles);

% --- Executes on button press in GroundMultipathCheck.
function GroundMultipathCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundMultipathCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundMultipathCheck
UpdateView(handles);

% --- Executes on button press in GroundShadowCheck.
function GroundShadowCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundShadowCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundShadowCheck
UpdateView(handles);

% --- Executes on button press in GroundNorthCheck.
function GroundNorthCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundNorthCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundNorthCheck
UpdateView(handles);

% --- Executes on button press in TurnAllOff.
%function TurnAllOff_Callback(hObject, eventdata, handles)
% hObject    handle to TurnAllOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%UpdateView(handles);

% --- Executes on button press in GroundBorderCheck.
function GroundBorderCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundBorderCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundBorderCheck
UpdateView(handles);

% --- Executes on button press in SlantNorthCheck.
function SlantNorthCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlantNorthCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlantNorthCheck
UpdateView(handles);

% --- Executes on button press in SlantShadowCheck.
function SlantShadowCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlantShadowCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlantShadowCheck
UpdateView(handles);

% --- Executes on button press in SlantMultipathCheck.
function SlantMultipathCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlantMultipathCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlantMultipathCheck
UpdateView(handles);

% --- Executes on button press in SlantLayoverCheck.
function SlantLayoverCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlantLayoverCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlantLayoverCheck
UpdateView(handles);

% --- Executes on button press in SlantNormalCheck.
function SlantNormalCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SlantNormalCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlantNormalCheck
UpdateView(handles);

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

%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

%get iq filename
[fname, pathstr] = uigetfile( sar_file_extensions( 'complex' ),...
    'Open Image File',pathstr,'MultiSelect', 'off');
if isnumeric(fname)
    return;
end

filename = [pathstr fname];
reader_obj = open_reader(filename);
if iscell(reader_obj)
    reader_obj = reader_obj{1};
end
handles.meta = reader_obj.get_meta();

try
    MetaIcon_Complex(handles.meta,'handle',handles.metaicon);
catch
end

UpdateView(handles);

setpref('matlab_sar_toolbox','last_used_directory',pathstr); %store path

function UpdateView(handles)

%get settings from GUI
PlotOptions.VelocityCheck = get(handles.VelocityCheck,'Value');
PlotOptions.RangeCheck = get(handles.RangeCheck,'Value');
PlotOptions.RangeGroundCheck = get(handles.RangeGroundCheck,'Value');
PlotOptions.SlantPlaneCheck = get(handles.SlantPlaneCheck,'Value');
PlotOptions.GroundPlaneCheck = get(handles.GroundPlaneCheck,'Value');
PlotOptions.GroundBorderCheck = get(handles.GroundBorderCheck,'Value');
PlotOptions.GroundSlantCheck = get(handles.GroundSlantCheck,'Value');
PlotOptions.SSPCheck = get(handles.SSPCheck,'Value');
PlotOptions.GroundNormalCheck = get(handles.GroundNormalCheck,'Value');
PlotOptions.GroundLayoverCheck = get(handles.GroundLayoverCheck,'Value');
PlotOptions.GroundMultipathCheck = get(handles.GroundMultipathCheck,'Value');
PlotOptions.GroundShadowCheck = get(handles.GroundShadowCheck,'Value');
PlotOptions.GroundNorthCheck = get(handles.GroundNorthCheck,'Value');
PlotOptions.SlantNormalCheck = get(handles.SlantNormalCheck,'Value');
PlotOptions.SlantLayoverCheck = get(handles.SlantLayoverCheck,'Value');
PlotOptions.SlantMultipathCheck = get(handles.SlantMultipathCheck,'Value');
PlotOptions.SlantShadowCheck = get(handles.SlantShadowCheck,'Value');
PlotOptions.SlantNorthCheck = get(handles.SlantNorthCheck,'Value');
PlotOptions.GrazeCheck = get(handles.Graze,'Value');
PlotOptions.AzimuthCheck = get(handles.Azimuth,'Value');
PlotOptions.SlopeCheck = get(handles.Slope,'Value');
PlotOptions.TwistCheck = get(handles.Twist,'Value');
PlotOptions.SquintGCheck = get(handles.SquintG,'Value');
PlotOptions.SquintSCheck = get(handles.SquintS,'Value');
PlotOptions.DCACheck = get(handles.DCACheck,'Value');
    
%delete anything in the image and legend axis
kids = get(handles.image,'Children');
delete(kids);

%this gets rid of anything in the flex legend...the number below may need
%to be modified as more controls are added to the main figure.  Just
%deleting the handle that is passed back from the flexlegend is
%insufficient as if the user clicks too fast it can't keep up and there are
%orphan legend entries...
try
    kids = get(handles.figure1,'Children');
    for i=1:length(kids)-9
        kids2 = get(kids(i),'Children');
        delete(kids2);
    end
catch
end

handles.Vectors = Plot3DImageGeometry(handles.meta,PlotOptions,handles.image);

% Update handles structure
guidata(handles.figure1, handles);

% --- Executes on button press in SSPCheck.
function SSPCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SSPCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SSPCheck
UpdateView(handles);

% --- Executes on button press in SlantOff.
function SlantOff_Callback(hObject, eventdata, handles)
% hObject    handle to SlantOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(handles.SlantOff,'String'),'All Off')
    set(handles.SlantNormalCheck,'Value',0);
    set(handles.SlantLayoverCheck,'Value',0);
    set(handles.SlantMultipathCheck,'Value',0);
    set(handles.SlantShadowCheck,'Value',0);
    set(handles.SlantNorthCheck,'Value',0);
    set(handles.SlantOff,'String','All On');
else
    set(handles.SlantNormalCheck,'Value',1);
    set(handles.SlantLayoverCheck,'Value',1);
    set(handles.SlantMultipathCheck,'Value',1);
    set(handles.SlantShadowCheck,'Value',1);
    set(handles.SlantNorthCheck,'Value',1);
    set(handles.SlantOff,'String','All Off');
end

UpdateView(handles);

% --- Executes on button press in GroundOff.
function GroundOff_Callback(hObject, eventdata, handles)
% hObject    handle to GroundOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(handles.GroundOff,'String'),'All Off')
    set(handles.GroundNormalCheck,'Value',0);
    set(handles.GroundLayoverCheck,'Value',0);
    set(handles.GroundMultipathCheck,'Value',0);
    set(handles.GroundShadowCheck,'Value',0);
    set(handles.GroundNorthCheck,'Value',0);
    set(handles.GroundOff,'String','All On');
else
    set(handles.GroundNormalCheck,'Value',1);
    set(handles.GroundLayoverCheck,'Value',1);
    set(handles.GroundMultipathCheck,'Value',1);
    set(handles.GroundShadowCheck,'Value',1);
    set(handles.GroundNorthCheck,'Value',1);
    set(handles.GroundOff,'String','All Off');
end

UpdateView(handles);


% --- Executes on selection change in PerspectiveCombo.
function PerspectiveCombo_Callback(hObject, eventdata, handles)
% hObject    handle to PerspectiveCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PerspectiveCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PerspectiveCombo
contents = cellstr(get(hObject,'String'));
View = contents{get(hObject,'Value')};

CamPos = get(handles.image,'CameraPosition');
if strcmpi(View,'Slant Plane Normal')
    CamPos = handles.Vectors.SPN*norm(CamPos);
    set(handles.image,'CameraTargetMode','auto');
    set(handles.image,'CameraUpVectorMode','auto');
elseif strcmpi(View,'Ground Plane Normal')
    CamPos = handles.Vectors.GPN*norm(CamPos);
    set(handles.image,'CameraTarget',[0 0 0]);
    if handles.meta.SCPCOA.SideOfTrack == 'R'
        set(handles.image,'CameraUpVector',[0 1 0]);
    else
        set(handles.image,'CameraUpVector',[0 -1 0]);
    end
elseif strcmpi(View,'Range Vector')
    CamPos = handles.Vectors.R*norm(CamPos);
    set(handles.image,'CameraTargetMode','auto');
    set(handles.image,'CameraUpVectorMode','auto');
end
set(handles.image,'CameraPosition',CamPos);


% --- Executes during object creation, after setting all properties.
function PerspectiveCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PerspectiveCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function ApPercent_Callback(hObject, eventdata, handles)
% hObject    handle to ApPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApPercent as text
%        str2double(get(hObject,'String')) returns contents of ApPercent as a double


% --- Executes during object creation, after setting all properties.
function ApPercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApPercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrentTime_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentTime as text
%        str2double(get(hObject,'String')) returns contents of CurrentTime as a double


% --- Executes during object creation, after setting all properties.
function CurrentTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Graze.
function Graze_Callback(hObject, eventdata, handles)
% hObject    handle to Graze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Graze
UpdateView(handles);

% --- Executes on button press in Azimuth.
function Azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Azimuth
UpdateView(handles);

% --- Executes on button press in Slope.
function Slope_Callback(hObject, eventdata, handles)
% hObject    handle to Slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Slope
UpdateView(handles);

% --- Executes on button press in Twist.
function Twist_Callback(hObject, eventdata, handles)
% hObject    handle to Twist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Twist
UpdateView(handles);


% --- Executes on button press in SquintG.
function SquintG_Callback(hObject, eventdata, handles)
% hObject    handle to SquintG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SquintG
UpdateView(handles);

% --- Executes on button press in SquintS.
function SquintS_Callback(hObject, eventdata, handles)
% hObject    handle to SquintS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SquintS
UpdateView(handles);


% --- Executes on button press in DCACheck.
function DCACheck_Callback(hObject, eventdata, handles)
% hObject    handle to DCACheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DCACheck
UpdateView(handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

typ = get( gcf, 'SelectionType' );    
if strcmp(typ,'alt')
    %use Preset View on Right Button Click
    contents = cellstr(get(handles.PerspectiveCombo,'String'));
    View = contents{get(handles.PerspectiveCombo,'Value')};
    CamPos = get(handles.image,'CameraPosition');
    if strcmpi(View,'Slant Plane Normal')
        CamPos = handles.Vectors.SPN*norm(CamPos);
        set(handles.image,'CameraTargetMode','auto');
        set(handles.image,'CameraUpVectorMode','auto');
    elseif strcmpi(View,'Ground Plane Normal')
        CamPos = handles.Vectors.GPN*norm(CamPos);
        set(handles.image,'CameraTarget',[0 0 0]);
        set(handles.image,'CameraUpVector',[0 1 0]);
    elseif strcmpi(View,'Range Vector')
        CamPos = handles.Vectors.R*norm(CamPos);
        set(handles.image,'CameraTargetMode','auto');
        set(handles.image,'CameraUpVectorMode','auto');
    end
    set(handles.image,'CameraPosition',CamPos);
  return;
end

handles.buttondown = 1;
handles.camerapos = get(handles.image,'CameraPosition');
set(handles.image,'CameraTargetMode','auto');
set(handles.image,'CameraUpVectorMode','auto');
handles.pos = get(gcf,'CurrentPoint');

guidata(hObject, handles);


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.ChangingView == 1)
    return;
end

if (handles.buttondown == 1)    
    handles.ChangingView = 1;
    guidata(hObject, handles);
    
    pos = get(gcf,'CurrentPoint');
    delta = handles.pos-pos;
       
    CamPos = handles.camerapos;
    CamPosN = norm(CamPos);
    Elevation = asind(CamPos(3)/CamPosN);
    Az = atan2d(CamPos(2),CamPos(1));
    Az = Az + delta(1)*2;
    Elevation = Elevation + delta(2)*2;
    if Elevation > 90
        Elevation = 90-(Elevation-90);
        %Az = Az+180;
    elseif Elevation < -90
        Elevation = -90-(Elevation+90);
        %Az = Az+180;
    end
            
    CamPos(3) = sind(Elevation);
    CamPos(1) = cosd(Az)*(1-abs(CamPos(3))^2);
    CamPos(2) = sind(Az)*(1-abs(CamPos(3))^2);
    
    CamPos = CamPos./norm(CamPos);
    CamPos = CamPos*CamPosN;
    set(handles.image,'CameraPosition',CamPos);
    
    handles.ChangingView = 0;
    guidata(hObject, handles);
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.buttondown = 0;
% Update handles structure
guidata(hObject, handles);
