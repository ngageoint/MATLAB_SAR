function varargout = TaserPreferences(varargin)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TaserPreferences_OpeningFcn, ...
                   'gui_OutputFcn',  @TaserPreferences_OutputFcn, ...
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


% --- Executes just before TaserPreferences is made visible.
function TaserPreferences_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TaserPreferences (see VARARGIN)

% Choose default command line output for TaserPreferences
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

Prefs = varargin{1};
if Prefs.ImageOverview
    set(handles.DisplayImageOverview,'Value',1);
else
    set(handles.DisplayMPSD,'Value',1);
end

if strcmp(Prefs.SICDLocation,'Colocated')
    set(handles.Colocated,'Value',1);
elseif strcmp(Prefs.SICDLocation,'TempFolder')
    set(handles.TempFolder,'Value',1);
else
    set(handles.SICDFolder,'Value',1);
    set(handles.SICDFolderDir,'String',Prefs.SICDLocation);
end
set(handles.MaxRes,'String',Prefs.MaxResolution);
set(handles.UseOverview,'Value',Prefs.UseOverview);
set(handles.NumPulses,'String',Prefs.NumPulses);
set(handles.ProcessingPulses,'String',Prefs.ProcessingPulses);
set(handles.dBCheck,'Value',Prefs.dBCheck);
set(handles.ComplexDefaultCheck,'Value',Prefs.ComplexDefaultCheck);
set(handles.PHDDefaultCheck,'Value',Prefs.PHDDefaultCheck);
set(handles.ComplexTemplate,'String',Prefs.ComplexTemplate);
set(handles.PHDTemplate,'String',Prefs.PHDTemplate);

handles.InitPrefs = Prefs;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TaserPreferences wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = TaserPreferences_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Prefs;
delete(handles.figure1);

% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Prefs.ImageOverview = get(handles.DisplayImageOverview,'Value');
if get(handles.Colocated,'Value')
    handles.Prefs.SICDLocation = 'Colocated';
elseif get(handles.TempFolder,'Value')
    handles.Prefs.SICDLocation = 'TempFolder';
else
    handles.Prefs.SICDLocation = get(handles.SICDFolderDir,'String');
end
handles.Prefs.MaxResolution = str2double(get(handles.MaxRes,'String'));
handles.Prefs.UseOverview = get(handles.UseOverview,'Value');
handles.Prefs.NumPulses = str2double(get(handles.NumPulses,'String'));
handles.Prefs.ProcessingPulses = str2double(get(handles.ProcessingPulses,'String'));
handles.Prefs.dBCheck = get(handles.dBCheck,'Value');
handles.Prefs.ComplexDefaultCheck = get(handles.ComplexDefaultCheck,'Value');
handles.Prefs.PHDDefaultCheck = get(handles.PHDDefaultCheck,'Value');
handles.Prefs.ComplexTemplate = get(handles.ComplexTemplate,'String');
handles.Prefs.PHDTemplate = get(handles.PHDTemplate,'String');

% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1);

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Prefs = handles.InitPrefs;
% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1);

function NumPulses_Callback(hObject, eventdata, handles)
% hObject    handle to NumPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumPulses as text
%        str2double(get(hObject,'String')) returns contents of NumPulses as a double


% --- Executes during object creation, after setting all properties.
function NumPulses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dBCheck.
function dBCheck_Callback(hObject, eventdata, handles)
% hObject    handle to dBCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dBCheck



function ProcessingPulses_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessingPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ProcessingPulses as text
%        str2double(get(hObject,'String')) returns contents of ProcessingPulses as a double


% --- Executes during object creation, after setting all properties.
function ProcessingPulses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProcessingPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


% --- Executes on button press in PHDDefaultCheck.
function PHDDefaultCheck_Callback(hObject, eventdata, handles)
% hObject    handle to PHDDefaultCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PHDDefaultCheck



function PHDTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to PHDTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PHDTemplate as text
%        str2double(get(hObject,'String')) returns contents of PHDTemplate as a double


% --- Executes during object creation, after setting all properties.
function PHDTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PHDTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BrowsePHDTemplate.
function BrowsePHDTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to BrowsePHDTemplate (see GCBO)
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

[FileName,PathName] = uigetfile(pathstr);
fname = [PathName filesep FileName];
set(handles.PHDTemplate,'String',fname);

% --- Executes on button press in ComplexDefaultCheck.
function ComplexDefaultCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ComplexDefaultCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ComplexDefaultCheck



function ComplexTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to ComplexTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ComplexTemplate as text
%        str2double(get(hObject,'String')) returns contents of ComplexTemplate as a double


% --- Executes during object creation, after setting all properties.
function ComplexTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ComplexTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseComplexTemplate.
function BrowseComplexTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseComplexTemplate (see GCBO)
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

[FileName,PathName] = uigetfile(pathstr);
fname = [PathName filesep FileName];
set(handles.ComplexTemplate,'String',fname);



function MaxRes_Callback(hObject, eventdata, handles)
% hObject    handle to MaxRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxRes as text
%        str2double(get(hObject,'String')) returns contents of MaxRes as a double


% --- Executes during object creation, after setting all properties.
function MaxRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SICDFolderDir_Callback(hObject, eventdata, handles)
% hObject    handle to SICDFolderDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SICDFolderDir as text
%        str2double(get(hObject,'String')) returns contents of SICDFolderDir as a double


% --- Executes during object creation, after setting all properties.
function SICDFolderDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SICDFolderDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseDir.
function BrowseDir_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseDir (see GCBO)
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

folder_name = uigetdir(pathstr);
set(handles.SICDFolderDir,'String',folder_name);


% --- Executes on button press in UseOverview.
function UseOverview_Callback(hObject, eventdata, handles)
% hObject    handle to UseOverview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UseOverview
