function varargout = SARLS(varargin)
% SARLS MATLAB code for SARLS.fig
%      SARLS, by itself, creates a new SARLS or raises the existing
%      singleton*.
%
%      H = SARLS returns the handle to a new SARLS or the handle to
%      the existing singleton*.
%
%      SARLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SARLS.M with the given input arguments.
%
%      SARLS('Property','Value',...) creates a new SARLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SARLS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SARLS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SARLS

% Last Modified by GUIDE v2.5 17-Jul-2020 08:54:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SARLS_OpeningFcn, ...
                   'gui_OutputFcn',  @SARLS_OutputFcn, ...
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


% --- Executes just before SARLS is made visible.
function SARLS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SARLS (see VARARGIN)

% Choose default command line output for SARLS
handles.output = hObject;

%read file extensions file and populate table
DirName = fileparts(which('SARLS'));
ExtSettings = jsondecode(fileread(fullfile(DirName, 'FileExtensions.json')));
set(handles.ExtTable,'ColumnName',{'Extension','Toggle'});
set(handles.ExtTable,'ColumnWidth',{90 60});
set(handles.ExtTable,'ColumnFormat',{'char' 'logical'});
set(handles.ExtTable,'ColumnEditable',[false true]);
data = fieldnames(ExtSettings);
data(:,2) = num2cell(structfun(@(x) x, ExtSettings));
set(handles.ExtTable,'Data',data);

%read meta fields file and populate table headings
MetaFields = jsondecode(fileread(fullfile(DirName, 'MetaFields.json')));
handles.fields = fieldnames(MetaFields);
handles.FieldFlags = structfun(@(x) double(x), MetaFields);
set(handles.MetaTable,'ColumnName',handles.fields(logical(handles.FieldFlags)));
set(handles.MetaTable,'Data',[]);

set(handles.RecursiveCheck,'Value',1);
set(handles.Process,'Enable','off');
set(handles.NumImagesText,'Visible','off');
set(handles.CSVFile,'Enable','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SARLS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SARLS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function DirName_Callback(hObject, eventdata, handles)
% hObject    handle to DirName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirName as text
%        str2double(get(hObject,'String')) returns contents of DirName as a double


% --- Executes during object creation, after setting all properties.
function DirName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirName (see GCBO)
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

if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

DirName = uigetdir(pathstr);

set(handles.DirName,'String',DirName);
set(handles.Process,'Enable','on');

handles.filenames = GetFileNames(handles);

% Update handles structure
guidata(hObject, handles);

function files = GetFileNames(handles)

DirName = get(handles.DirName,'String');

%get list of files based on selected extensions and recursive check
data = get(handles.ExtTable,'data');
exts = data([data{:,2}],1);
if ~isempty(exts) && ~isempty(DirName)
    if get(handles.RecursiveCheck,'Value')
        files = rdir(DirName,exts);
    else
        %get files in current directory with matching extensions   
        temp = cellfun(@(x) dir(fullfile(DirName, ['*.' x])), exts, 'UniformOutput', false);
        files = arrayfun(@(x) fullfile(DirName, x.name), [temp{:}], 'UniformOutput', false);
    end
else 
    files = {};
end

set(handles.NumImagesText,'Visible','on');
set(handles.NumImagesText,'String',[num2str(length(files)) ' Potential SAR files']);

% --- Executes on button press in RecursiveCheck.
function RecursiveCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RecursiveCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.filenames = GetFileNames(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Process.
function Process_Callback(hObject, eventdata, handles)
% hObject    handle to Process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get data for all fields, then just display selected fields, that way if
%fields are added, the metadata doesn;t need to be recomputed (potentially time consuming) 
handles.metadata = GetMetaData(handles.filenames,handles.fields);

% Update handles structure
guidata(hObject, handles);

DisplayMetaTable(handles);

set(handles.CSVFile,'Enable','on');

% --- Executes on button press in TableFields.
function TableFields_Callback(hObject, eventdata, handles)
% hObject    handle to TableFields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.fields,handles.FieldFlags] = SelectMetaFields(handles.fields,handles.FieldFlags);

% Update handles structure
guidata(hObject, handles);

if isfield(handles,'metadata')
    DisplayMetaTable(handles);
end

function DisplayMetaTable(handles)
%DISPLAYMETATABLE Updates MetaTable data

column_names = handles.fields(logical(handles.FieldFlags));
data = cell(size(handles.metadata,2), numel(column_names));
for i = 1:numel(column_names)
    data(:,i) = {handles.metadata.(column_names{i})};
end

set(handles.MetaTable,'ColumnName', column_names);
set(handles.MetaTable,'data', data);

% --- Executes when entered data in editable cell(s) in ExtTable.
function ExtTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ExtTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

handles.filenames = GetFileNames(handles);

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CSVFile_Callback(hObject, eventdata, handles)
% hObject    handle to CSVFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

%get CSV filename
[fname, pathstr] = uiputfile('*.csv','Save CSV File',pathstr);

%get column headers and data from MetaTable
Headers = get(handles.MetaTable,'ColumnName');
data = get(handles.MetaTable,'data');

M = vertcat(Headers',data);

cell2csv([pathstr filesep fname],M); 

