function varargout = Image2KMLGUI(varargin)
% IMAGE2KMLGUI M-file for Image2KMLGUI.fig
%      IMAGE2KMLGUI, by itself, creates a new IMAGE2KMLGUI or raises the existing
%      singleton*.
%
%      H = IMAGE2KMLGUI returns the handle to a new IMAGE2KMLGUI or the handle to
%      the existing singleton*.
%
%      IMAGE2KMLGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE2KMLGUI.M with the given input arguments.
%
%      IMAGE2KMLGUI('Property','Value',...) creates a new IMAGE2KMLGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Image2KMLGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Image2KMLGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Edit the above text to modify the response to help Image2KMLGUI

% Last Modified by GUIDE v2.5 22-Jul-2014 15:37:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Image2KMLGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Image2KMLGUI_OutputFcn, ...
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


% --- Executes just before Image2KMLGUI is made visible.
function Image2KMLGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Image2KMLGUI (see VARARGIN)

% Choose default command line output for Image2KMLGUI
handles.output = hObject;

set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');
set(handles.metaicon,'Color',[51/255 102/255 153/255]);

set(handles.GroundOverlayCheck,'Value',1);
handles.DecType = 'none';
handles.DecSize = 3000;

set(handles.PlotSRP,'Value',1);
set(handles.PlotPath,'Value',1);
set(handles.PlotGroundTrack,'Value',1);
set(handles.PlotCollectionWedge,'Value',1);
set(handles.PlotCenterApVector,'Value',1);
set(handles.PlotPRFLines,'Value',1);

% Update handles structure
guidata(hObject, handles);

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.parse(varargin{:});

if ~isempty(p.Results.filename)
    handles = LoadImage(p.Results.filename,handles);
end


% --- Outputs from this function are returned to the command line.
function varargout = Image2KMLGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseFile.
function BrowseFile_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseFile (see GCBO)
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
[fname, path] = uigetfile(sar_file_extensions({'complex','phd'}),'Open SAR Data File',pathstr);
if isnumeric(fname)
    return;
end

setpref('matlab_sar_toolbox','last_used_directory',path); %store path

filename = strcat(path,fname);
handles = LoadImage(filename,handles);

% Update handles structure
guidata(hObject, handles);


function handles = LoadImage(filename,handles)

set(handles.filename,'String',filename);

% determine if input file is complex image or phase history data
if ~isempty(guess_ph_format(filename))
    MetaIcon_PHD(filename,'handle',handles.metaicon);
    % Should disable all of the image controls here (overlay, boundary, and fill)
else
    MetaIcon_Complex(filename,'handle',handles.metaicon);
end


% --- Executes on button press in WriteKML.
function WriteKML_Callback(hObject, eventdata, handles)
% hObject    handle to WriteKML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get kml name from user
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

%open input file
filename = get(handles.filename,'String');
if ~exist(filename,'file')
    return;
end
is_phase_history = ~isempty(guess_ph_format(filename));
if is_phase_history
    reader_obj = open_ph_reader(filename);
else
    reader_obj = open_reader(filename);
end
if ~iscell(reader_obj) % Assure cell array format
    reader_obj = {reader_obj};
end
meta = reader_obj{1}.get_meta();
[path_ignore, fnameNoExt] = fileparts(filename);
% Open output file
if get(handles.GroundOverlayCheck,'Value') 
    [fname, path] = uiputfile( {'*.kmz','KMZ Files (*.kmz)'},'Save KMZ File',[pathstr fnameNoExt]);
else
    [fname, path] = uiputfile( {'*.kml','KML Files (*.kml)'},'Save KML File',[pathstr fnameNoExt]);
end
if ~fname, return; end; % Cancel was pressed
if ~isfield(meta, 'CollectionInfo') && isfield(meta, 'CollectionID')
    meta.CollectionInfo = meta.CollectionID;  % CPHD
end
k = kml(meta.CollectionInfo.CoreName); %create an kmltoolbox object
k.filename = fullfile(path, fname);
k.zip = get(handles.GroundOverlayCheck,'Value'); % Only zip if we have overlay images

%set border values
if get(handles.PlotImageBorder,'Value')
    BorderColor = 'Red';
    FillColor = 'Red';
else
    BorderColor = [];
    FillColor = [];
end

%loop through segments
NumSegments = length(reader_obj);  
if NumSegments>1
    wb = waitbar(0);
end
for i=1:NumSegments
    if NumSegments>1
        waitbar(i/NumSegments,wb,sprintf('Processing segment %d of %d',i,NumSegments));
    end
    
    if get(handles.GroundOverlayCheck,'Value')
        %image filename is state with path relative to KML
        ImName = sprintf('Overlay_%d.png',i); % Will be put in same directory where KML is
    else
        ImName = '';
    end
    add_sar_2kml(k,reader_obj{i},...
        'srp',get(handles.PlotSRP,'Value'),...
        'border_color',BorderColor,...
        'border_thickness',1,...
        'fill_color',FillColor,...
        'transparency',0.25,...
        'overlay_max_size',handles.DecSize,...
        'overlay_decimate',handles.DecType,...
        'overlay_filename',ImName,...
        'sensor_path',get(handles.PlotPath,'Value'),...
        'sensor_ground_track',get(handles.PlotGroundTrack,'Value'),...
        'collection_wedge',get(handles.PlotCollectionWedge,'Value'),...
        'center_aperture_vector',get(handles.PlotCenterApVector,'Value'),...
        'ambiguity_bounds',get(handles.PlotPRFLines,'Value'),...
        'description_function',@(x) default_SICD_description(x,get(handles.English,'Value')));
end

k.save();
if ~isempty(k.includeFiles)
    delete(k.includeFiles{:});
end
reader_obj{1}.close; % Multiple segments, but only a single file to close
if NumSegments>1
    close(wb); % Waitbar
end

msgbox('Done creating shapefile');


% --- Executes on button press in GroundOverlayCheck.
function GroundOverlayCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundOverlayCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GroundOverlayCheck

if get(handles.GroundOverlayCheck,'Value')
    %set size of output overlay image based on Quality specified
    switch get(handles.OverlayQualityCombo,'Value');
        case 1
            handles.DecSize = 5000;
        case 2
            handles.DecSize = 3000;
        otherwise
            handles.DecSize = 1000;
    end           

    %remap, max decimated typically looks better with "brighter" remap
    if get(handles.MaxDecimateCheck,'Value')
        handles.DecType = 'max';
    else
        handles.DecType = 'none';
    end
else
    handles.DecSize = 0;
    handles.DecType = '';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function OverlayQualityCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverlayQualityCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotSRP.
function PlotSRP_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSRP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotSRP


% --- Executes on button press in PlotPath.
function PlotPath_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotPath


% --- Executes on button press in PlotGroundTrack.
function PlotGroundTrack_Callback(hObject, eventdata, handles)
% hObject    handle to PlotGroundTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotGroundTrack


% --- Executes on button press in PlotCollectionWedge.
function PlotCollectionWedge_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCollectionWedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotCollectionWedge


% --- Executes on button press in PlotCenterApVector.
function PlotCenterApVector_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCenterApVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotCenterApVector


% --- Executes on button press in PlotPRFLines.
function PlotPRFLines_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPRFLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotPRFLines


% --- Executes on button press in Batch.
function Batch_Callback(hObject, eventdata, handles)
% hObject    handle to Batch (see GCBO)
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

data_dir = uigetdir(pathstr, 'Select directory of SAR data');
if ~data_dir, return; end; % Cancel
setpref('matlab_sar_toolbox','last_used_directory',data_dir); %store path

if get(handles.GroundOverlayCheck,'Value') 
    [fname, path] = uiputfile( {'*.kmz','KMZ Files (*.kmz)'},'Save KMZ File',pathstr);
else
    [fname, path] = uiputfile( {'*.kml','KML Files (*.kml)'},'Save KML File',pathstr);
end
if ~fname, return; end; % Cancel was pressed
[path_ignore, fname_no_extension] = fileparts(fname);

%set border values
if get(handles.PlotImageBorder,'Value')
    BorderColor = 'Red';
    FillColor = 'Red';
else
    BorderColor = [];
    FillColor = [];
end

image2kmlbatch(data_dir, strcat(path,fname), ...
    'name',fname_no_extension,...
    'srp',get(handles.PlotSRP,'Value'),...
    'border_color',BorderColor,...
    'border_thickness',1,...
    'fill_color',FillColor,...
    'transparency',0.25,...
    'overlay_max_size',handles.DecSize,...
    'overlay_decimate',handles.DecType,...
    'sensor_path',get(handles.PlotPath,'Value'),...
    'sensor_ground_track',get(handles.PlotGroundTrack,'Value'),...
    'collection_wedge',get(handles.PlotCollectionWedge,'Value'),...
    'center_aperture_vector',get(handles.PlotCenterApVector,'Value'),...
    'ambiguity_bounds',get(handles.PlotPRFLines,'Value'),...
    'description_function',@(x) default_SICD_description(x,get(handles.English,'Value')));
msgbox('Done creating shapefile');


% --- Executes on selection change in OverlayQualityCombo.
function OverlayQualityCombo_Callback(hObject, eventdata, handles)
% hObject    handle to OverlayQualityCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OverlayQualityCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OverlayQualityCombo


% --- Executes on button press in PlotImageBorder.
function PlotImageBorder_Callback(hObject, eventdata, handles)
% hObject    handle to PlotImageBorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotImageBorder

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////