function varargout = SelectMetaFields(varargin)
% SELECTMETAFIELDS MATLAB code for SelectMetaFields.fig
%      SELECTMETAFIELDS, by itself, creates a new SELECTMETAFIELDS or raises the existing
%      singleton*.
%
%      H = SELECTMETAFIELDS returns the handle to a new SELECTMETAFIELDS or the handle to
%      the existing singleton*.
%
%      SELECTMETAFIELDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTMETAFIELDS.M with the given input arguments.
%
%      SELECTMETAFIELDS('Property','Value',...) creates a new SELECTMETAFIELDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectMetaFields_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectMetaFields_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectMetaFields

% Last Modified by GUIDE v2.5 15-Jul-2020 14:54:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectMetaFields_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectMetaFields_OutputFcn, ...
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


% --- Executes just before SelectMetaFields is made visible.
function SelectMetaFields_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectMetaFields (see VARARGIN)

% Choose default command line output for SelectMetaFields
handles.output = hObject;

handles.fields = varargin{1};
handles.FieldFlags = varargin{2};

%populate available and selected fields
Count1 = 0;
Count2 = 0;
AvailFields = {};
SelFields = {};
for i=1:length(handles.fields)
    if handles.FieldFlags(i)
        Count2 = Count2 + 1;
        SelFields{Count2} = handles.fields{i};
    else
        Count1 = Count1 + 1;
        AvailFields{Count1} = handles.fields{i};
    end
end

set(handles.AvailableFields,'String',AvailFields);
set(handles.SelectedFields,'String',SelFields);

if isempty(AvailFields)
    set(handles.AddField,'Enable','off');
    set(handles.AddAll,'Enable','off');
end

if isempty(SelFields)
    set(handles.RemoveField,'Enable','off');
    set(handles.RemoveAll,'Enable','off');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectMetaFields wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SelectMetaFields_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.fields;
varargout{2} = handles.FieldFlags;
delete(handles.figure1);

% --- Executes on selection change in AvailableFields.
function AvailableFields_Callback(hObject, eventdata, handles)
% hObject    handle to AvailableFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AvailableFields contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AvailableFields


% --- Executes during object creation, after setting all properties.
function AvailableFields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AvailableFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SelectedFields.
function SelectedFields_Callback(hObject, eventdata, handles)
% hObject    handle to SelectedFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectedFields contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectedFields


% --- Executes during object creation, after setting all properties.
function SelectedFields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectedFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddField.
function AddField_Callback(hObject, eventdata, handles)
% hObject    handle to AddField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get selection, remove from available fields and add to selected fields
Strings = get(handles.AvailableFields,'String');
Val = get(handles.AvailableFields,'Value');
StringsNew = {};
count = 0;
for i=1:length(Strings)
    if i~=Val
        count = count + 1;
        StringsNew{count} = Strings{i};
    end
end
set(handles.AvailableFields,'Value',1);
set(handles.AvailableFields,'String',StringsNew);
if isempty(StringsNew)
    set(handles.AddField,'Enable','off');
    set(handles.AddAll,'Enable','off');
else
    set(handles.AddField,'Enable','on');
    set(handles.AddAll,'Enable','on');
end

set(handles.RemoveField,'Enable','on');
set(handles.RemoveAll,'Enable','on');

Strings2 = get(handles.SelectedFields,'String');
Strings2{end+1} = Strings{Val};
set(handles.SelectedFields,'String',Strings2);

% --- Executes on button press in RemoveField.
function RemoveField_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get selection, remove from available fields and add to selected fields
Strings = get(handles.SelectedFields,'String');
Val = get(handles.SelectedFields,'Value');
StringsNew = {};
count = 0;
for i=1:length(Strings)
    if i~=Val
        count = count + 1;
        StringsNew{count} = Strings{i};
    end
end
set(handles.SelectedFields,'Value',1);
set(handles.SelectedFields,'String',StringsNew);
if isempty(StringsNew)
    set(handles.RemoveField,'Enable','off');
    set(handles.RemoveAll,'Enable','off');
else
    set(handles.RemoveField,'Enable','on');
    set(handles.RemoveAll,'Enable','on');
end
set(handles.AddField,'Enable','on');
set(handles.AddAll,'Enable','on');

Strings2 = get(handles.AvailableFields,'String');
Strings2{end+1} = Strings{Val};
set(handles.AvailableFields,'String',Strings2);

% --- Executes on button press in AddAll.
function AddAll_Callback(hObject, eventdata, handles)
% hObject    handle to AddAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Strings = get(handles.AvailableFields,'String');
Strings2 = get(handles.SelectedFields,'String');

set(handles.AvailableFields,'Value',1);
set(handles.AvailableFields,'String',[]);
set(handles.SelectedFields,'String',vertcat(Strings2,Strings));

set(handles.AddField,'Enable','off');
set(handles.AddAll,'Enable','off');

set(handles.RemoveField,'Enable','on');
set(handles.RemoveAll,'Enable','on');

% --- Executes on button press in RemoveAll.
function RemoveAll_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Strings = get(handles.AvailableFields,'String');
Strings2 = get(handles.SelectedFields,'String');

set(handles.SelectedFields,'Value',1);
set(handles.SelectedFields,'String',[]);
set(handles.AvailableFields,'String',vertcat(Strings,Strings2));

set(handles.RemoveField,'Enable','off');
set(handles.RemoveAll,'Enable','off');

set(handles.AddField,'Enable','on');
set(handles.AddAll,'Enable','on');

% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update fields and FieldFlags based on the list box contents
SelFields = get(handles.SelectedFields,'String');
AvailFields = get(handles.AvailableFields,'String');

handles.fields = vertcat(SelFields,AvailFields);
handles.FieldFlags = [ones(length(SelFields),1); zeros(length(AvailFields),1)];

guidata(hObject, handles);

uiresume(handles.figure1);

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1);
