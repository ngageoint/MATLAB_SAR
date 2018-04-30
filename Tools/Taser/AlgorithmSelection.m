function varargout = AlgorithmSelection(varargin)
% ALGORITHMSELECTION MATLAB code for AlgorithmSelection.fig
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
%      ALGORITHMSELECTION, by itself, creates a new ALGORITHMSELECTION or raises the existing
%      singleton*.
%
%      H = ALGORITHMSELECTION returns the handle to a new ALGORITHMSELECTION or the handle to
%      the existing singleton*.
%
%      ALGORITHMSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALGORITHMSELECTION.M with the given input arguments.
%
%      ALGORITHMSELECTION('Property','Value',...) creates a new ALGORITHMSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AlgorithmSelection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AlgorithmSelection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AlgorithmSelection

% Last Modified by GUIDE v2.5 13-Jun-2011 07:36:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AlgorithmSelection_OpeningFcn, ...
                   'gui_OutputFcn',  @AlgorithmSelection_OutputFcn, ...
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


% --- Executes just before AlgorithmSelection is made visible.
function AlgorithmSelection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AlgorithmSelection (see VARARGIN)

% Choose default command line output for AlgorithmSelection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

Algs = varargin{1};

SetAlgorithmLists(handles,Algs);

handles.Algorithms = Algs;
handles.InitialAlgs = Algs;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AlgorithmSelection wait for user response (see UIRESUME)
uiwait(handles.figure1);


function SetAlgorithmLists(handles,Algorithms)

SelectedCount = 0;
AvailableCount = 0;
for i=1:length(Algorithms)
    if strcmp(Algorithms(i).Selected,'Yes')
        SelectedCount = SelectedCount + 1;
        SelectedNames{SelectedCount} = Algorithms(i).TextName;
    else
        AvailableCount = AvailableCount + 1;
        AvailableNames{AvailableCount} = Algorithms(i).TextName;
    end
end

if (AvailableCount > 0)
    set(handles.AvailableAlgorithms,'Value',1);
    set(handles.AvailableAlgorithms,'String',AvailableNames);
else
    set(handles.AvailableAlgorithms,'Value',0);
    set(handles.AvailableAlgorithms,'String',[]);
end

if (SelectedCount > 0)
    set(handles.SelectedAlgorithms,'Value',1);
    set(handles.SelectedAlgorithms,'String',SelectedNames);
else
    set(handles.SelectedAlgorithms,'Value',0);
    set(handles.SelectedAlgorithms,'String',[]);
end

if (SelectedCount > 0)
    set(handles.MoveLeft,'enable','on');
else
    set(handles.MoveLeft,'enable','off');
end
if (AvailableCount > 0)
    set(handles.MoveRight,'enable','on');
else
    set(handles.MoveRight,'enable','off');
end

% --- Outputs from this function are returned to the command line.
function varargout = AlgorithmSelection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Algorithms;
delete(handles.figure1);

% --- Executes on selection change in AvailableAlgorithms.
function AvailableAlgorithms_Callback(hObject, eventdata, handles)
% hObject    handle to AvailableAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AvailableAlgorithms contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AvailableAlgorithms


% --- Executes during object creation, after setting all properties.
function AvailableAlgorithms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AvailableAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SelectedAlgorithms.
function SelectedAlgorithms_Callback(hObject, eventdata, handles)
% hObject    handle to SelectedAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectedAlgorithms contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectedAlgorithms


% --- Executes during object creation, after setting all properties.
function SelectedAlgorithms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectedAlgorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MoveRight.
function MoveRight_Callback(hObject, eventdata, handles)
% hObject    handle to MoveRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get selection from Available List
Index = get(handles.AvailableAlgorithms,'Value');
Strings = get(handles.AvailableAlgorithms,'String');
Algorithm = Strings{Index};

%find Algorithm and change status
for i=1:length(handles.Algorithms)
    if strcmp(Algorithm,handles.Algorithms(i).TextName)
        handles.Algorithms(i).Selected = 'Yes';
    end
end

SetAlgorithmLists(handles,handles.Algorithms)

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in MoveLeft.
function MoveLeft_Callback(hObject, eventdata, handles)
% hObject    handle to MoveLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get selection from Available List
Index = get(handles.SelectedAlgorithms,'Value');
Strings = get(handles.SelectedAlgorithms,'String');
Algorithm = Strings{Index};

%find Algorithm and change status
for i=1:length(handles.Algorithms)
    if strcmp(Algorithm,handles.Algorithms(i).TextName)
        handles.Algorithms(i).Selected = 'No';
    end
end

SetAlgorithmLists(handles,handles.Algorithms)

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Algorithms = handles.InitialAlgs;
guidata(hObject, handles);
uiresume(handles.figure1);

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////