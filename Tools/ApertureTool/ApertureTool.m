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
%     - initial version, based on previous focus tools written by Tim Cox.
%       MITM tool and readers written by Wade Schwartzkopf.  Initial
%       concept based on tool written by Matt Banta.
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
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
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

% Setup image and PHD axes
remap_list = getremaplist();
default_remap = find(strcmpi(remap_list,'densityremap'),1,'first');
set(handles.remap_combo,'String',remap_list);
set(handles.remap_combo,'Value',default_remap);
image_rep_list = {'amplitude' 'Pauli' 'AlphaEntropy'};
set(handles.image_rep_combo,'String',image_rep_list);
set(handles.image_rep_combo,'Value',1);
old_img_position = getpixelposition(handles.image);
handles.imghandle = imshow(0, 'Parent', handles.image);
handles.hSP1 = imscrollpanel(handles.figure1,handles.imghandle); 
setpixelposition(handles.hSP1,old_img_position);
handles.apiSP1 = iptgetapi(handles.hSP1);

%main tab
set(handles.AlwaysApply,'Value',1);

%animation tab
set(handles.SlowTime,'Value',1);
set(handles.NumFrames,'String',7);
set(handles.ApertureFraction,'String',0.25);
set(handles.FrameRate,'String',5);
set(handles.CycleCheck,'Value',1);
set(handles.AlphaFactor,'String',0.20);
set(handles.Percentile,'String',50);
set(handles.CMapCutoff,'String',150);

%disable controls until image is loaded
set([handles.FirstFrame,...
    handles.PrevFrame,...
    handles.Play,...
    handles.Pause,...
    handles.NextFrame,...
    handles.LastFrame],'enable','off');
set(handles.image_rep_combo,'value',1,'visible','off');
set(handles.polarization_uipanel,'visible','off');

%resolution tab
set(handles.Resolution,'Value',1);
set(handles.MaxMin,'Value',1);
%disable controls until image is loaded
set([handles.MinRes,...
    handles.MaxRes,...
    handles.StepSize,...
    handles.FirstFrameRes,...
    handles.PrevFrameRes,...
    handles.PlayRes,...
    handles.PauseRes,...
    handles.NextFrameRes,...
    handles.LastFrameRes],'enable','off');

%save image tab
set([handles.SaveImage handles.SaveMetaIcon],'enable','off');

%save movie tab
set(handles.MovieFrameRate,'String',5);
set(handles.ModeAnnotate,'Value',1);
set(handles.FrameAnnotate,'Value',1);
set(handles.ResAnnotate,'Value',1);
set(handles.AvgResAnnotate,'Value',0);
%disable controls until image is loaded
set(handles.Record,'enable','off');
set(handles.StopRecord,'enable','off');

%add radio-button callbacks
set(handles.ResModeRadio,'SelectionChangeFcn',@ResModeChange);
set(handles.DomainPanel,'SelectionChangeFcn',@DomainChange);
set(handles.FilterPanel,'SelectionChangeFcn',@FilterChange);

%set status flags
handles.ResMode = 0;
handles.Recording = 0;

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
    LoadImage(p.Results.filename,hObject,handles,AOI,segment);
end

% --- Outputs from this function are returned to the command line.
function varargout = ApertureTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ImageContrast.
function ImageContrast_Callback(hObject, eventdata, handles)
% hObject    handle to ImageContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imcontrast(handles.image);


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
updatePHD_axesLabels(handles)


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


% --- Executes on selection change in ImageZoomCombo.
function ImageZoomCombo_Callback(hObject, eventdata, handles)
% hObject    handle to ImageZoomCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


function CrossRangeStart_Callback(hObject, eventdata, handles)
% hObject    handle to CrossRangeStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossRangeStart as text
%        str2double(get(hObject,'String')) returns contents of CrossRangeStart as a double


% --- Executes during object creation, after setting all properties.
function CrossRangeStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossRangeStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RangeStart_Callback(hObject, eventdata, handles)
% hObject    handle to RangeStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RangeStart as text
%        str2double(get(hObject,'String')) returns contents of RangeStart as a double


% --- Executes during object creation, after setting all properties.
function RangeStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RangeStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CrossRangeStop_Callback(hObject, eventdata, handles)
% hObject    handle to CrossRangeStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossRangeStop as text
%        str2double(get(hObject,'String')) returns contents of CrossRangeStop as a double


% --- Executes during object creation, after setting all properties.
function CrossRangeStop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossRangeStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RangeStop_Callback(hObject, eventdata, handles)
% hObject    handle to RangeStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RangeStop as text
%        str2double(get(hObject,'String')) returns contents of RangeStop as a double


% --- Executes during object creation, after setting all properties.
function RangeStop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RangeStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CrossRangeFraction_Callback(hObject, eventdata, handles)
% hObject    handle to CrossRangeFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossRangeFraction as text
%        str2double(get(hObject,'String')) returns contents of CrossRangeFraction as a double


% --- Executes during object creation, after setting all properties.
function CrossRangeFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossRangeFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RangeFraction_Callback(hObject, eventdata, handles)
% hObject    handle to RangeFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RangeFraction as text
%        str2double(get(hObject,'String')) returns contents of RangeFraction as a double


% --- Executes during object creation, after setting all properties.
function RangeFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RangeFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CrossRangeRes_Callback(hObject, eventdata, handles)
% hObject    handle to CrossRangeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossRangeRes as text
%        str2double(get(hObject,'String')) returns contents of CrossRangeRes as a double


% --- Executes during object creation, after setting all properties.
function CrossRangeRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossRangeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RangeRes_Callback(hObject, eventdata, handles)
% hObject    handle to RangeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RangeRes as text
%        str2double(get(hObject,'String')) returns contents of RangeRes as a double


% --- Executes during object creation, after setting all properties.
function RangeRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RangeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CrossRangeSS_Callback(hObject, eventdata, handles)
% hObject    handle to CrossRangeSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossRangeSS as text
%        str2double(get(hObject,'String')) returns contents of CrossRangeSS as a double


% --- Executes during object creation, after setting all properties.
function CrossRangeSS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossRangeSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RangeSS_Callback(hObject, eventdata, handles)
% hObject    handle to RangeSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RangeSS as text
%        str2double(get(hObject,'String')) returns contents of RangeSS as a double


% --- Executes during object creation, after setting all properties.
function RangeSS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RangeSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ImageName_Callback(hObject, eventdata, handles)
% hObject    handle to ImageName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageName as text
%        str2double(get(hObject,'String')) returns contents of ImageName as a double


% --- Executes during object creation, after setting all properties.
function ImageName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseImage.
function BrowseImage_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


LoadImage('',hObject,handles);



% --- Executes on button press in UniformWeightingCheck.
function UniformWeightCheck_Callback(hObject, eventdata, handles)
% hObject    handle to UniformWeightingCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UniformWeightingCheck
process_phd(hObject,handles);


% --- Executes on button press in CenterPHDCheckAz.
function CenterPHDCheckAz_Callback(hObject, eventdata, handles)
% hObject    handle to CenterPHDCheckAz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CenterPHDCheckAz
if get(hObject,'Value') % In general, we can't deskew in range and cross-range simultaneously
    if ~isscalar(handles.meta.Grid.Col.DeltaKCOAPoly) || ...
            ~isscalar(handles.meta.Grid.Row.DeltaKCOAPoly)
        set(handles.CenterPHDCheckRng,'Value',false);
    end
    handles.manual_offset = [0 0];
elseif all(handles.meta.Grid.Row.DeltaKCOAPoly(:)==0)&&handles.isvalid_row_deltakcoapoly
    set(handles.CenterPHDCheckRng,'Value',true);
end
process_phd(hObject,handles);


% --- Executes on button press in CenterPHDCheckRng.
function CenterPHDCheckRng_Callback(hObject, eventdata, handles)
% hObject    handle to CenterPHDCheckRng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CenterPHDCheckRng
if get(hObject,'Value') % In general, we can't deskew in range and cross-range simultaneously
    if ~isscalar(handles.meta.Grid.Col.DeltaKCOAPoly) || ...
            ~isscalar(handles.meta.Grid.Row.DeltaKCOAPoly)
        set(handles.CenterPHDCheckAz,'Value',false);
    end
    handles.manual_offset = [0 0];
elseif all(handles.meta.Grid.Col.DeltaKCOAPoly(:)==0) && handles.isvalid_col_deltakcoapoly
    set(handles.CenterPHDCheckAz,'Value',true);
end
process_phd(hObject,handles);


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


function MaxLoadSize_Callback(hObject, eventdata, handles)
% hObject    handle to MaxLoadSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxLoadSize as text
%        str2double(get(hObject,'String')) returns contents of MaxLoadSize as a double


% --- Executes during object creation, after setting all properties.
function MaxLoadSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxLoadSize (see GCBO)
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
newpos(pos, hObject);


% Load AOI from file
function LoadImage(filenames,hObject,handles,varargin)

if ((nargin<1)||isempty(filenames)) % If no filename was give, use dialog box to ask for one
    % Recall last interactively selected path used
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
    if(iscell(filenames)), filenames=sort(filenames); end;
    setpref('matlab_sar_toolbox','last_used_directory',pathname);
    if isfield(handles,'phdhandle')
        handles = rmfield(handles,'phdhandle');
    end
else % Path was already passed in with filenames
    pathname='';
end
fullfilename = {}; % Reset if files were previously open
if(iscell(filenames)) % Multiple files requested
    for j=1:length(filenames)
        fullfilename{j}=[pathname filenames{j}];
    end
elseif(filenames)
    fullfilename{1}=[pathname filenames];
else % filename=0.  Cancel was pressed, instead of a file being chosen.
    return;
end
if ~isempty(pathname)
    setpref('matlab_sar_toolbox','last_used_directory',pathname); %store path
end

% Open readers for all files, and build cell array of all readers
reader_obj ={};
sf_ind = []; fn_ind = []; % Indices to filenames and subimages for each reader object
for i = 1:length(fullfilename)
    newreader = open_reader(fullfilename{i});
    fn_ind = [fn_ind i*ones(size(newreader))];
    sf_ind = [sf_ind 1:numel(newreader)];
    reader_obj = [reader_obj newreader];
end
clear newreader;

set(handles.ImageName,'String',fullfilename{1});
% reader_obj=open_reader(fullfilename);

[readerobj_indices,reader_obj]=group_reader_objs_by_pol(reader_obj); % Consolidate reader objects into polarimetric sets

%get metadata
if (iscell(reader_obj))
    meta=reader_obj{1}.get_meta(); % Assume first frame for getting default frame/aoi size
    % Adjust later if frame changes
else
    meta=reader_obj.get_meta();
end
nx = double(meta.ImageData.NumCols);
ny = double(meta.ImageData.NumRows);

% Select AOI and frame number if necessary
if isempty(varargin) || isempty(varargin{1}) || varargin{1}(1) == 0
    %determine if we'll just load the image or if we'll load an AOI
    %selection MITM_Viewer
    MaxLoadSize = str2double(get(handles.MaxLoadSize,'String'));
    if (nx <= MaxLoadSize && ny <= MaxLoadSize)
        %load image directly into tool
        aoi = [1 1 nx ny];
        FrameNum = 1;
    else
        %launch a MITM viewer with an AOI selection box
         ButtonName = questdlg(sprintf('%s %i x %i\n%s %i x %i\n%s',...
            'Suggested AOI dimensions are less than',MaxLoadSize,MaxLoadSize,...
            'Current AOI dimensions are',nx,ny,...
            'Do you want to redefine the AOI?'), ...
            'AptertureTool AOI', ...
            'Yes', 'No', 'Yes');
        if isempty(ButtonName), return; end; %window was closed CANCEL
        if strcmpi(ButtonName,'No')
            aoi = [1 1 nx ny]; %use AOI anyway
            FrameNum = 1;
        else
            [aoi, FrameNum] = mitm_viewer(fullfilename,'mode','aoi','closeAfterSelect',true);
            FrameNum = find(cellfun(@(x) isequal(FrameNum,sf_ind(x)), readerobj_indices));
        end
        if isempty(aoi), return; end; % Cancel was selected
    end
else
    aoi = varargin{1};
    if nargin>4 && ~isempty(varargin{2}) && all(varargin{2} ~= 0)
        FrameNum = find(cellfun(@(x) isequal(varargin{2},sf_ind(x)), readerobj_indices));
    else
        FrameNum = 1;
    end
end
if (iscell(reader_obj))
    meta=reader_obj{FrameNum}.get_meta(); % Assume first frame for getting datasize
end
    
try % If some metadata for icon is not given, don't make tool fail
    MetaIcon_Complex(meta,'handle',handles.metaicon);
    set(handles.SaveMetaIcon,'enable','on');
catch
    cla(handles.metaicon);
end

% Handle complex data normalization parameters from metadata
% Set data normalization controls based on metadata
if isfield(meta,'Grid') && all(isfield(meta.Grid,{'Row','Col'})) && ...
        isfield(meta.Grid.Row,'WgtType') && ...
        isfield(meta.Grid.Col,'WgtType') && ...
        ((strcmpi(meta.Grid.Col.WgtType,'UNIFORM') && ...
        strcmpi(meta.Grid.Row.WgtType,'UNIFORM')) || ...
        (isfield(meta.Grid.Row.WgtType,'WindowName') && ...
        isfield(meta.Grid.Col.WgtType,'WindowName') && ...
        strcmpi(meta.Grid.Col.WgtType.WindowName,'UNIFORM') && ...
        strcmpi(meta.Grid.Row.WgtType.WindowName,'UNIFORM')))
    % If we already know that data is uniform weighted, this is not an option.
    set(handles.UniformWeightingCheck,'Value',true);
    set(handles.UniformWeightingCheck,'Enable','off');
else
    set(handles.UniformWeightingCheck,'Value',false);
    set(handles.UniformWeightingCheck,'Enable','on');
end
% Set default values for required frequency support deskew values
handles.isvalid_row_deltakcoapoly = isfield(meta,'Grid') && ...
    isfield(meta.Grid,'Row') && isfield(meta.Grid.Row,'DeltaKCOAPoly') && ...
    ~isempty(meta.Grid.Row.DeltaKCOAPoly) && ...
    all(isfinite(meta.Grid.Row.DeltaKCOAPoly(:)));
if handles.isvalid_row_deltakcoapoly && any(meta.Grid.Row.DeltaKCOAPoly(:)~=0)
    set(handles.CenterPHDCheckRng,'Enable','on');
else
    set(handles.CenterPHDCheckRng,'Enable','off');
    meta.Grid.Row.DeltaKCOAPoly = 0; % Give default value (no shift)
end
handles.isvalid_col_deltakcoapoly = isfield(meta,'Grid') && ...
    isfield(meta.Grid,'Col') && isfield(meta.Grid.Col,'DeltaKCOAPoly') && ...
    ~isempty(meta.Grid.Col.DeltaKCOAPoly) && ...
    all(isfinite(meta.Grid.Col.DeltaKCOAPoly(:)));
if handles.isvalid_col_deltakcoapoly && any(meta.Grid.Col.DeltaKCOAPoly(:)~=0)
    set(handles.CenterPHDCheckAz,'Enable','on');
else
    set(handles.CenterPHDCheckAz,'Enable','off');
    meta.Grid.Col.DeltaKCOAPoly = 0; % Give default value (no shift)
end
% Is it possible to deskew frequency support in both dimensions at once?
if isscalar(meta.Grid.Col.DeltaKCOAPoly) && ...
        isscalar(meta.Grid.Row.DeltaKCOAPoly)
    set(handles.CenterPHDCheckAz,'Value',handles.isvalid_col_deltakcoapoly); % Center both dimensions
    set(handles.CenterPHDCheckRng,'Value',handles.isvalid_row_deltakcoapoly);
else % Make user select dimension to deskew
    set(handles.CenterPHDCheckAz,'Value',handles.isvalid_col_deltakcoapoly && ...
        all(meta.Grid.Col.DeltaKCOAPoly(:) == 0));
    set(handles.CenterPHDCheckRng,'Value',handles.isvalid_row_deltakcoapoly && ...
        all(meta.Grid.Row.DeltaKCOAPoly(:) == 0));
end
handles.manual_offset = [0 0];
% Determine weighting function
try
    handles.weight_fun_az = sicdweight2fun(meta.Grid.Col);
    handles.weight_fun_rng = sicdweight2fun(meta.Grid.Row);
    handles.weighting_computed_from_data = false;
catch % Unable to determine weighting functions from metadata.
    % Estimate weighting from complex data instead
    if iscell(reader_obj)
        handles.weight_fun_az = estimate_weighting_file(reader_obj{FrameNum}, 1);
        handles.weight_fun_rng = estimate_weighting_file(reader_obj{FrameNum}, 2);
    else
        handles.weight_fun_az = estimate_weighting_file(reader_obj, 1);
        handles.weight_fun_rng = estimate_weighting_file(reader_obj, 2);
    end
    handles.weighting_computed_from_data = true;
end
% Determine if we have metadata to compute zeropad accurately
if isfield(meta,'Grid') && all(isfield(meta.Grid,{'Row','Col'})) && ...
        all(isfield(meta.Grid.Row,{'SS','ImpRespBW'})) && ...
        all(isfield(meta.Grid.Col,{'SS','ImpRespBW'}))
    handles.AzPad = max(1,1/(meta.Grid.Col.SS*meta.Grid.Col.ImpRespBW));
    handles.RnPad = max(1,1/(meta.Grid.Row.SS*meta.Grid.Row.ImpRespBW));
else % Otherwise, just default to full data.
    handles.AzPad = 1;
    handles.RnPad = 1;
end % User can always manual adjust zeropad area

% Use correct FFT function based on data
% Sgn of -1 means that a forward FFT is the correct transform to move from
% image domain to spatial frequency.  However, we want to flip the
% direction, so that high RF frequencies are at the top (low indices) when
% displayed in MATLAB, and time/angle increases toward the right (higher
% indices).  So we swap the FFT/IFFTs which, in effect, flips left/right
% and up/down in spatial frequency.
if isfield(meta,'Grid') && all(isfield(meta.Grid,{'Row','Col'})) && ...
        isfield(meta.Grid.Row,'Sgn') && isfield(meta.Grid.Col,'Sgn') && ...
        meta.Grid.Row.Sgn > 0 && meta.Grid.Col.Sgn > 0
    handles.fft_sp = @fft2; % Function to transform from image to spatial frequency
    handles.fft_im = @ifft2; % Function to transfrom from spatial frequency to image
elseif isfield(meta,'Grid') && all(isfield(meta.Grid,{'Row','Col'})) && ...
        isfield(meta.Grid.Row,'Sgn') && isfield(meta.Grid.Col,'Sgn') && ...
        meta.Grid.Row.Sgn ~= meta.Grid.Col.Sgn
    error('APERTURETOOL:MIXED_FFT_SGN','Mixed FFT signs not handled.');
else
    handles.fft_sp = @ifft2; % Function to transform from image to spatial frequency
    handles.fft_im = @fft2; % Function to transfrom from spatial frequency to image
end

%get polarization information
if isfield(meta,'ImageFormation') && isfield(meta.ImageFormation,'TxRcvPolarizationProc')
    handles.TxRcvPolarizationProc = meta.ImageFormation.TxRcvPolarizationProc;
    if ~iscell(handles.TxRcvPolarizationProc)
        handles.TxRcvPolarizationProc = {handles.TxRcvPolarizationProc};
        set(handles.image_rep_combo,'value',1)
    else
        set(handles.image_rep_combo,'value',...
            find(strcmpi('Pauli',get(handles.image_rep_combo,'string'))))
    end
    handles.polar_PHD_idx = 1;
    if numel(handles.TxRcvPolarizationProc) == 1
        set(handles.image_rep_combo,'value',1,'visible','off')
        set(handles.polarization_uipanel,'visible','off')
    else
        set(handles.image_rep_combo,'value',2,'visible','on')
        set(handles.polarization_uipanel,'visible','on')
    end
else
    set(handles.image_rep_combo,'value',1,'visible','off')
    set(handles.polarization_uipanel,'visible','off')
end
updatePolarToggle(handles)

% Read data
handles.complex_data = double(reader_obj{FrameNum}.read_chip(...
    [aoi(1) aoi(1)+aoi(3)-1],[aoi(2) aoi(2)+aoi(4)-1]));
reader_obj{FrameNum}.close();
% Condition data so that imcontrast control works better
handles.complex_data = 10*handles.complex_data./mean(handles.complex_data(:));
handles.aoi = aoi; % store aoi for later
handles.segment = FrameNum;

%set chip size in Load File Tab
set(handles.nx,'String',aoi(3));
set(handles.ny,'String',aoi(4));

% This is derived redundant information, but we compute it once, so we
% don't need to do this in multiple places:
handles.ZpLimsAz = round((aoi(3)/2) * (1 - 1/handles.AzPad)) + 1;
handles.ZpLimsAz(2) = aoi(3) - handles.ZpLimsAz + 1;
handles.ZpWidthAz = diff(handles.ZpLimsAz)+1;
handles.ZpLimsRn = round((aoi(4)/2) * (1 - 1/handles.RnPad)) + 1;
handles.ZpLimsRn(2) = aoi(4) - handles.ZpLimsRn + 1;
handles.ZpWidthRn = diff(handles.ZpLimsRn)+1;

% update metadata
% chip info
meta.ImageData.FirstRow = meta.ImageData.FirstRow + aoi(2) - 1;
meta.ImageData.FirstCol = meta.ImageData.FirstCol + aoi(1) - 1;
meta.ImageData.NumRows = aoi(4);
meta.ImageData.NumCols = aoi(3);
try
    meta = add_sicd_corners(meta);
end
% If ImpRespWid is not defined, but ImpRespBW is, we just use what it would
% be for uniform weighting (since we're probably going to want to apply a
% uniform weighting to the data anyway.)
set(handles.PHDZoomCombo,'value',1); % Reset to "Full"
resolution_defined = true;
if ~isfield(meta.Grid.Col,'ImpRespWid')
    if isfield(meta.Grid.Col,'ImpRespBW')
        meta.Grid.Col.ImpRespWid = .886/meta.Grid.Col.ImpRespBW;
    else % Just set to zero if we don't know it
        meta.Grid.Col.ImpRespWid = 0;
        meta.Grid.Col.SS = 0;
        resolution_defined = false;
    end
end
if ~isfield(meta.Grid.Row,'ImpRespWid')
    if isfield(meta.Grid.Row,'ImpRespBW')
        meta.Grid.Row.ImpRespWid = .886/meta.Grid.Row.ImpRespBW;
    else % Just set to zero if we don't know it
        meta.Grid.Row.ImpRespWid = 0;
        meta.Grid.Row.SS = 0;
        resolution_defined = false;
    end
end
if resolution_defined
    set(handles.AperturePercent,'enable','on');
    set(handles.GroundProjCheck,'enable','on');
    set(handles.Resolution,'enable','on');
else % Force aperture percent resolution mode
    set(handles.AperturePercent,'Value',1);
    set(handles.AperturePercent,'enable','off');
    set(handles.GroundProjCheck,'Value',0);
    set(handles.GroundProjCheck,'enable','off');
    set(handles.Resolution,'enable','off');
    set(handles.MinResUnits,'String','Percent');
    set(handles.MaxResUnits,'String','Percent');
    set(handles.StepSizeUnits,'String','Percent');
end
% Without twist and graze, we can't calculate ground resolution.
if ~isfield(meta,'SCPCOA') || ~isfield(meta.SCPCOA,'TwistAng') ||...
    ~isfield(meta.SCPCOA,'GrazeAng')
    set(handles.GroundProjCheck,'enable','off');
    set(handles.GroundProjCheck,'value',false);

    meta.SCPCOA.TwistAng = 0;  % Resort to slant plane resolution.
    meta.SCPCOA.GrazeAng = 0;  % Resort to slant plane resolution.
    warning('ApertureTool:NoGroundResolution',...
        ['Insufficient metadata to compute ground resolution.  ',...
        'Slant plane resolution will be used if available.']);
end
handles.meta = meta; % store metadata for later

process_phd(hObject,handles); % Apply PHD options
handles = guidata(hObject); % process_phd change handles (We may fix this later.)

%set initial magnifications to fit
mag = handles.apiSP1.findFitMag();
handles.apiSP1.setMagnification(mag);
set(handles.ImageZoomCombo,'Value',1);
strings = get(handles.ImageZoomCombo,'String');
strings{1} = sprintf('Fit (%d%%)',round(mag*100));
set(handles.ImageZoomCombo,'String',strings);


% Draw imrect on phd
if isfield(handles,'H') && ishandle(handles.H)
    delete(handles.H); % Delete any previous imrects
end
handles.H = imrect(handles.phd, [handles.ZpLimsAz(1), handles.ZpLimsRn(1),...
    handles.ZpWidthAz handles.ZpWidthRn]);
if verLessThan('matlab', '7.6')
    disp('Warning. Matlab version does not support imrect constraint function');
else
    setPositionConstraintFcn(handles.H,...
        makeConstrainToRectFcn('imrect',handles.ZpLimsAz,handles.ZpLimsRn));
end
api_handle = iptgetapi(handles.H);
api_handle.addNewPositionCallback(@(x) newpos(x,hObject));

% Now enable relevant controls
set(handles.Play,'enable','on');
set(handles.LastFrame,'enable','on');
set(handles.FirstFrame,'enable','on');
set(handles.NextFrame,'enable','off'); % Full aperture is not a frame, so
set(handles.PrevFrame,'enable','off'); % no next or previous
set(handles.SaveImage,'enable','on');
set(handles.Record,'enable','on');
set(handles.ManualAdjust,'enable','on');

% Set up default settings for resolution tab
SetResAnimationSettings(handles);

% Sample spacing only needs to be set at load time, since it never changes
is_english = get(handles.UnitsCheck,'Value');
set_res_fields(meta.Grid.Col.SS, handles.CrossRangeSS, handles.CRSSUnits, is_english);
set_res_fields(meta.Grid.Row.SS, handles.RangeSS, handles.SSUnits, is_english);

% Update handles structure
guidata(hObject, handles);

% Update resolution fields
newpos(getPosition(handles.H),hObject);


% Condition PHD as specified (center, deweight, left/right flip)
function process_phd(hObject,handles)
for ii = 1:size(handles.complex_data,3) %treat PHD for each polarization separately
    cdata = handles.complex_data(:,:,ii);
    meta = handles.meta;
    if all(handles.manual_offset==0)
        % If possible deweight in an already deskewed dimension before deskewing in
        % the other dimension.  We will deweight the other dimension after deskew.
        if all(meta.Grid.Col.DeltaKCOAPoly(:) == 0) && ...
                all(meta.Grid.Row.DeltaKCOAPoly(:) == 0) && ...
                get(handles.UniformWeightingCheck,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
            cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
        elseif all(meta.Grid.Col.DeltaKCOAPoly(:) == 0) && ...
                get(handles.UniformWeightingCheck,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
        elseif all(meta.Grid.Row.DeltaKCOAPoly(:) == 0) && ...
                get(handles.UniformWeightingCheck,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
        end
        % Frequency support deskew
        if get(handles.CenterPHDCheckAz,'Value')
            [ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( meta, 1 );
            if any(DeltaKCOAPoly(:)~=0) % Do we need to deskew frequency support?
                [cdata, new_DeltaKCOAPoly] = deskewmem(cdata, DeltaKCOAPoly, az_coords_m, rg_coords_m, 1, fft_sgn);
                % Deskewing might shift frequency support in other direction.
                if get(handles.UniformWeightingCheck,'Value')
                    cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
                end
            elseif ~get(handles.CenterPHDCheckRng,'Value') % If no deskew, at least attempt to recenter
                new_DeltaKCOAPoly = meta.Grid.Row.DeltaKCOAPoly;
            else
                new_DeltaKCOAPoly = 0;
            end
            if any(new_DeltaKCOAPoly(:)~=0)
                % In general, we cannot deskew in both directions at once.
                % However, we can at least recenter the non-deskewed dimension
                % with a simple shift (not a proper deskew, which is a spatially
                % variant shift).
                deltaKCOA = sicd_polyval2d(new_DeltaKCOAPoly,...
                    az_coords_m(round(numel(az_coords_m)/2)),...
                    rg_coords_m(round(numel(rg_coords_m)/2))); % Get shift at center of AOI
                % deltaKCOA is scalar so it results in a spatially invariant
                % shift which will not affect the deskew diresction.
                cdata = deskewmem(cdata, deltaKCOA, az_coords_m, rg_coords_m, 2, meta.Grid.Row.Sgn);
            end
        end
        if get(handles.CenterPHDCheckRng,'Value')
            [ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( meta, 2 );
            if any(DeltaKCOAPoly(:)~=0) % Do we need to center frequency support?
                [cdata, new_DeltaKCOAPoly] = deskewmem(cdata, DeltaKCOAPoly, az_coords_m, rg_coords_m, 2, fft_sgn);
                if get(handles.UniformWeightingCheck,'Value')
                    cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
                end
            elseif ~get(handles.CenterPHDCheckAz,'Value') % If no deskew, at least attempt to recenter
                new_DeltaKCOAPoly = meta.Grid.Col.DeltaKCOAPoly;
            else
                new_DeltaKCOAPoly = 0;
            end
            if any(new_DeltaKCOAPoly(:)~=0)
                % See comments above for deskew and the resulting shift
                deltaKCOA = sicd_polyval2d(new_DeltaKCOAPoly,...
                    az_coords_m(round(numel(az_coords_m)/2)),...
                    rg_coords_m(round(numel(rg_coords_m)/2))); % Get shift at center of AOI
                cdata = deskewmem(cdata, deltaKCOA, az_coords_m, rg_coords_m, 1, meta.Grid.Col.Sgn);
            end
        end
    else
        cdata = handles.fft_im(circshift(handles.fft_sp(cdata),handles.manual_offset)); % Manual shift
        if get(handles.UniformWeightingCheck,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
            cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
        end
    end
    
    phd(:,:,ii) = fftshift(handles.fft_sp(cdata.'));
end
% Uniform weighting
if get(handles.UniformWeightingCheck,'Value') && ...
        isfield(meta.Grid.Row,'ImpRespBW') && isfield(meta.Grid.Col,'ImpRespBW')
    % Adjust impulse response width metadata for new weighting
    meta.Grid.Col.ImpRespWid = .886/meta.Grid.Col.ImpRespBW;
    meta.Grid.Row.ImpRespWid = .886/meta.Grid.Row.ImpRespBW;
end

handles.phasehistory = phd;

%store chip.  We will need to fix the contrast in resolution mode, so when
%it is selected, we will store the displayed sigmas above mean
% for ii = 1:size(phd,3)
%     chipmag(:,:,ii) = abs(handles.fft_im(phd(:,:,ii)));
% end
[chipmag,handles] = makeDisplayable(handles,phd);
handles.chip = chipmag;


%Draw Image 
remaps = get(handles.remap_combo,'String');
remap = remaps{get(handles.remap_combo,'Value')};
if strcmp(remap,'linearremap')
    chipmean = mean(chipmag(:));
    chipstd = std(chipmag(:));
    handles.apiSP1.replaceImage(chipmag,...
        'DisplayRange', [0 chipmean+3*chipstd], 'PreserveView', true);
else
    handles.apiSP1.replaceImage(chipmag, 'PreserveView', true);
end

%Draw PHD
if isfield(handles,'polar_PHD_idx')
    phdmag = abs(phd(:,:,handles.polar_PHD_idx));
else
    phdmag = abs(phd);
end
phdmean = mean(phdmag(:));
phdmag = 10.*phdmag./phdmean;
phdmean = 10;
phdstd = std(phdmag(:));

%if flight is left then the phd needs to be flipped horizontally to
%represent time going right (which we do in this country...)
if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SideOfTrack') && ...
        strcmp(meta.SCPCOA.SideOfTrack,'L')
    phdmag = fliplr(phdmag);
    flipXaxis = 0;
else
    flipXaxis = 1;
end

%make axes labels
try
    if isfield(handles,'phdYaxis')
        handles = rmfield(handles,'phdYaxis');
    end
    freq_width = (1/handles.meta.Grid.Row.SS)*(SPEED_OF_LIGHT/2);
    freq_ctr = handles.meta.Grid.Row.KCtr*(SPEED_OF_LIGHT/2);
    freq_limits = freq_ctr + ([-1 1]*freq_width/2);
    if isfield(handles.meta,'PFA') && isfield(handles.meta.PFA,'SpatialFreqSFPoly')
        freq_limits = freq_limits/handles.meta.PFA.SpatialFreqSFPoly(1);
    end
    if freq_ctr > 1e9
        freq_limits = freq_limits/1e9;
        handles.freqUnits = 'GHz';
    else
        freq_limits = freq_limits/1e6;
        handles.freqUnits = 'MHz';
    end
    handles.phdYaxis = linspace(freq_limits(2), freq_limits(1), size(phdmag,1));
end
try
    if isfield(handles,'phdXaxis')
        handles = rmfield(handles,'phdXaxis');
    end
    angle_width = (1/handles.meta.Grid.Col.SS) / handles.meta.Grid.Row.KCtr;
    if isfield(handles.meta.Grid.Col,'KCtr')
        angle_ctr = handles.meta.Grid.Col.KCtr;
    else
        angle_ctr = 0;
    end
    angle_limits = angle_ctr+([-1 1]*angle_width/2);
    if flipXaxis
        angle_limits = angle_limits([2 1]);
    end
    handles.phdXaxis = atand(linspace(angle_limits(1), angle_limits(2), size(phdmag,2)));
end

%update PHD image
if ~isfield(handles,'phdhandle')
   handles.phdhandle = imagesc(phdmag,'parent',handles.phd);
else
    set(handles.phdhandle,'cdata',phdmag)
end

handles.phdmagSize = size(phdmag);

updatePHD_axesLabels(handles)

% Update handles structure
guidata(hObject, handles);


function updatePHD_axesLabels(handles)
%update PHD axes labels

xlim = get(handles.phd,'xlim');
ylim = get(handles.phd,'ylim');

%setup axes labels
if isfield(handles,'phdXaxis')
    dxlim = diff(xlim);
    xticks = round(xlim(1))+[0 ceil(dxlim/4-1) ceil(dxlim/2-1) ...
        ceil(dxlim/4*3-1) dxlim];
    xticks(xticks < xlim(1)) = ceil(xlim(1));
    xticks(xticks > xlim(2)) = floor(xlim(2));
    xticklabels = cellfun(@(x) sprintf('%0.2g',x),num2cell(handles.phdXaxis(xticks)),'uniformoutput',false);
    xticklabels{3} = '0'; % Should be very close to zero, so this makes display cleaner.
    set(handles.phd,'xtick',xticks,'xticklabel',xticklabels);
    xlabel(handles.phd,'Polar Angle (degrees)','fontweight','bold')
else
    set(handles.phd, 'xtick', []);
end
if isfield(handles,'phdYaxis')
    dylim = diff(ylim);
    yticks = round(ylim(1))+[0 ceil(dylim/4-1) ceil(dylim/2-1) ...
        ceil(dylim/4*3-1) dylim];
    yticks(yticks < ylim(1)) = ceil(ylim(1));
    yticks(yticks > ylim(2)) = floor(ylim(2));
    set(handles.phd,'ytick',yticks,'yticklabel',...
        cellfun(@(x) sprintf('%0.4g',x),num2cell(handles.phdYaxis(yticks)),'uniformoutput',false));
    ylabel(handles.phd,sprintf('Frequency (%s)',handles.freqUnits),'fontweight','bold')
else
    set(handles.phd, 'ytick', []);
end


function [out,handles] = makeDisplayable(handles,phd)
% Do all processing to convert raw complex data into a real format
% MATLAB can display

for ii = 1:size(phd,3)
    out(:,:,ii) = handles.fft_im(phd(:,:,ii));
end

reps = get(handles.image_rep_combo,'String');
rep = reps{get(handles.image_rep_combo,'Value')};

if ~any(strcmpi(rep,{'amplitude' 'Pauli' 'AlphaEntropy'}))
    error('ApertureTool:ImageRep','ApertureTool: amplitude, Pauli, and AlphaEntropy are the only supported image representations');
end
if size(out,3) > 1
    if strcmpi(rep,'amplitude')
        %use single polarization channel (treat as single band)
        out = out(:,:,handles.polar_PHD_idx);
    else
        switch size(out,3)
            case 1 % Single band image; nothing to do
            case 2 % Dual-pol
                co_index=[find(strcmpi('H:H',handles.TxRcvPolarizationProc))...
                    find(strcmpi('V:V',handles.TxRcvPolarizationProc))];
                if length(co_index)>1 % Catch HH/VV
                    cross_index=co_index(2);
                    co_index=co_index(1);
                else
                    cross_index=[find(strcmpi('H:V',handles.TxRcvPolarizationProc))...
                        find(strcmpi('V:H',handles.TxRcvPolarizationProc))];
                end
                if strcmpi(rep,'AlphaEntropy')
                    [out(:,:,1),out(:,:,2),out(:,:,3)] =ComputeAlphaEntropy(out(:,:,co_index),...
                        out(:,:,cross_index),[],[]);
                else if strcmpi(rep,'Pauli')
                        out=cat(3,abs(out(:,:,co_index)),...
                            abs(out(:,:,cross_index)),...
                            abs(out(:,:,co_index)));
                    end
                end
            case 3 % RGB image; nothing to do
            case 4 % Quad-pol
                HH_index=find(strcmpi('H:H',handles.TxRcvPolarizationProc));
                HV_index=find(strcmpi('H:V',handles.TxRcvPolarizationProc));
                VH_index=find(strcmpi('V:H',handles.TxRcvPolarizationProc));
                VV_index=find(strcmpi('V:V',handles.TxRcvPolarizationProc));
                if strcmpi(rep,'AlphaEntropy')
                    [out(:,:,1), out(:,:,2), out(:,:,3)] =ComputeAlphaEntropy(out(:,:,HH_index),...
                        out(:,:,HV_index),out(:,:,VV_index),...
                        out(:,:,VH_index),1);
                    out(:,:,4) = [];
                else if strcmpi(rep,'Pauli')
                        out=cat(3,abs(out(:,:,HH_index)-out(:,:,VV_index)),...
                            abs(out(:,:,HV_index)+out(:,:,VH_index))/2,...
                            abs(out(:,:,HH_index)+out(:,:,VV_index)));
                    end
                end
        end
    end
else
    set(handles.image_rep_combo,'Value',find(strcmpi('amplitude',reps)));
end
handles.CurrentImage = out;

if isinteger(out) % Many remap functions won't work on complex ints
    out=single(out);
end

remaps = get(handles.remap_combo,'String');
remap = remaps{get(handles.remap_combo,'Value')};

% if nargin>2&&~isempty(remap)
% Apply remap per band
try
    for i=1:size(out,3)
        temp(:,:,i)=feval(remap,out(:,:,i));
    end
catch
    error('ApertureTool:InvalidRemap',['Error applying remap function: ' func2str(remap)]);
end
out=temp; % Use intermediate variable to allow for datatype change
if size(out,3)>1 && isfloat(out)% True color images must be between 0 and 1
    out=out-min(out(isfinite(out)));
    out=out/max(out(:));
end
% end


function newpos(pos, hObject)
handles = guidata(hObject);

%get filter type
if get(handles.FilterNone,'Value')
    WindowFilter = 0;
elseif get(handles.FilterGaussian,'Value')
    WindowFilter = 1;
elseif get(handles.Filterx4,'Value')
    WindowFilter = 2;
elseif get(handles.FilterHamming,'Value')
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

% Update GUI
% CRStart = round(((xmin-handles.ZpLimsAz(1)+1)/handles.ZpWidthAz)*100);
% CRStop = round(((xmax-handles.ZpLimsAz(1)+1)/handles.ZpWidthAz)*100);
% CRFrac = round(((xmax-xmin+1)/handles.ZpWidthAz)*100);
% RangeStart = round(((ymin-handles.ZpLimsRn(1)+1)/handles.ZpWidthRn)*100);
% RangeStop = round(((ymax-handles.ZpLimsRn(1)+1)/handles.ZpWidthRn)*100);
% RangeFrac = round(((ymax-ymin+1)/handles.ZpWidthRn)*100);
CRStart = ((xmin-handles.ZpLimsAz(1)+1)/handles.ZpWidthAz)*100;
CRStop = ((xmax-handles.ZpLimsAz(1)+1)/handles.ZpWidthAz)*100;
CRFrac = ((xmax-xmin+1)/handles.ZpWidthAz)*100;
RangeStart = ((ymin-handles.ZpLimsRn(1)+1)/handles.ZpWidthRn)*100;
RangeStop = ((ymax-handles.ZpLimsRn(1)+1)/handles.ZpWidthRn)*100;
RangeFrac = ((ymax-ymin+1)/handles.ZpWidthRn)*100;

%set up percentages/fractions etc
set(handles.CrossRangeStart,'String',CRStart);
set(handles.RangeStart,'String',RangeStart);
set(handles.CrossRangeStop,'String',CRStop);
set(handles.RangeStop,'String',RangeStop);
set(handles.CrossRangeFraction,'String',CRFrac);
set(handles.RangeFraction,'String',RangeFrac);

%compute resolutions
if get(handles.InverseSelection,'Value');
    RangeRes = handles.meta.Grid.Row.ImpRespWid./((100-RangeFrac)/100);
    AzRes = handles.meta.Grid.Col.ImpRespWid./((100-CRFrac)/100);
else
    RangeRes = handles.meta.Grid.Row.ImpRespWid./(RangeFrac/100);
    AzRes = handles.meta.Grid.Col.ImpRespWid./(CRFrac/100);
end
GroundRangeRes = RangeRes/cosd(handles.meta.SCPCOA.GrazeAng);
GroundAzRes = AzRes/cosd(handles.meta.SCPCOA.TwistAng);

is_english = get(handles.UnitsCheck,'Value');
set_res_fields(AzRes, handles.CrossRangeRes, handles.CRResUnits, is_english);
set_res_fields(RangeRes, handles.RangeRes, handles.ResUnits, is_english);
set_res_fields(GroundAzRes, handles.CrossRangeGroundRes, handles.CRGroundResUnits, is_english);
set_res_fields(GroundRangeRes, handles.RangeGroundRes, handles.GroundResUnits, is_english);

%need to modify xposition if flight is left
[ny, nx, nz] = size(handles.phasehistory);
if isfield(handles.meta,'SCPCOA') && isfield(handles.meta.SCPCOA,'SideOfTrack') && ...
        strcmp(handles.meta.SCPCOA.SideOfTrack,'L')
    xminold = xmin;
    xmin = round(nx-xmax);
    if (xmin<1);xmin=1;end;
    xmax = round(nx-xminold);    
end

%now reform image with selected aperture
if get(handles.InverseSelection,'Value')        
    Ap = handles.phasehistory;
    Ap(ymin:ymax,xmin:xmax,:) = 0;
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
    
    phdchip = FilterPHDChip(phdchip,WindowFilter,FilterRange,FilterAz);         
    
    Ap(ymin:ymax,xmin:xmax,:) = phdchip;    
end

[ImMag,handles] = makeDisplayable(handles,Ap);

remaps = get(handles.remap_combo,'String');
remap = remaps{get(handles.remap_combo,'Value')};
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

if (handles.Recording == 1)
    %record frame
    temp.Pos = pos;
    temp.CLim = get(handles.image,'CLim');
    temp.Mag = handles.apiSP1.getMagnification();
    R = handles.apiSP1.getVisibleImageRect();
    temp.Center = [round(R(1)+R(3)/2) round(R(2)+R(4)/2)]; 
    %save mode
    if isfield(handles, 'InLoop') && handles.InLoop
        %check for ResMode
        if get(handles.ResModeCheck,'Value')
            temp.Mode = 'Multi-Resolution';
        else
            %determine if this is slow of fast-time
            if get(handles.SlowTime,'Value')
                temp.Mode = 'Slow-Time';
            else
                temp.Mode = 'Fast-Time';
            end
        end
    else
        temp.Mode = 'User Interactive';
    end
        
    %save resolution (Ground Plane)
    temp.AzRes = str2double(get(handles.CrossRangeGroundRes,'String'));
    temp.RangeRes = str2double(get(handles.RangeGroundRes,'String'));
    temp.AzResUnits = get(handles.CRGroundResUnits,'String');
    temp.RangeResUnits = get(handles.GroundResUnits,'String');
    
    handles.aviparams = horzcat(handles.aviparams,temp);      
end

%if deployed, make sure image is updated now (do not want to do this
%normally, since it takes time)
if (isdeployed())
    drawnow expose;
end

% Update handles structure
guidata(handles.image, handles);

function phdchip = FilterPHDChip(phdchip,WindowFilter,FilterRange,FilterAz)   

[ny,nx] = size(phdchip);

switch WindowFilter        
    case 1 %gaussian
        if FilterRange
            win = gausswin(ny);
            phdchip = phdchip.*repmat(win,1,nx);
        end
        if FilterAz
            win = gausswin(nx);
            phdchip = phdchip.*repmat(win',ny,1);
        end
    case 2 %1/x^4
        if FilterRange
            win = x4win(ny);
            phdchip = phdchip.*repmat(win,1,nx);
        end
        if FilterAz
            win = x4win(nx);
            phdchip = phdchip.*repmat(win',ny,1);
        end
    case 3 %hamming
        if FilterRange
            win = hamming(ny);
            phdchip = phdchip.*repmat(win,1,nx);
        end
        if FilterAz
            win = hamming(nx);
            phdchip = phdchip.*repmat(win',ny,1);                
        end
    case 4 %cosine on pedastal
        if FilterRange
            win = cospedwin(ny,0.5);
            phdchip = phdchip.*repmat(win,1,nx);
        end
        if FilterAz
            win = cospedwin(nx,0.5);
            phdchip = phdchip.*repmat(win',ny,1);
        end
end


% Takes the appropiate MATLAB UIcontrols handles and sets them appropriately
function set_res_fields(value_meters, value_edit, units_text, is_english)
if is_english %english
    value = value_meters/FEET_TO_METERS;
    if (value < 1)
        value = value*12;
        set(units_text,'String','inches');
    else
        set(units_text,'String','feet');
    end
else %metric
    value = value_meters;
    if (value < 1)
        value = value*100;
        set(units_text,'String','cm');
    else
        set(units_text,'String','meters');
    end
end
set(value_edit,'String',sprintf('%.2f',value));


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


% --- Executes on button press in FirstFrame.
function FirstFrame_Callback(hObject, eventdata, handles)
% hObject    handle to FirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set([handles.NextFrame handles.LastFrame],'enable','on');
set([handles.FirstFrame handles.PrevFrame],'enable','off');
UpdateAnimationFrame(1, handles);


% --- Executes on button press in PrevFrame.
function PrevFrame_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%previous Frame
CurFrame = str2double(get(handles.CurFrame,'String')) - 1;

set([handles.NextFrame handles.LastFrame],'enable','on');

if (CurFrame == 1)    
    set([handles.FirstFrame handles.PrevFrame],'enable','off');
else
    set([handles.FirstFrame handles.PrevFrame],'enable','on');
end

UpdateAnimationFrame(CurFrame, handles);

% --- Executes on button press in Play.
function Play_Callback(hObject, eventdata, handles)
% hObject    handle to Play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get settings
NumFrames = str2double(get(handles.NumFrames,'String'));
FrameRate = str2double(get(handles.FrameRate,'String'));
CycleCheck = get(handles.CycleCheck,'Value'); % Cycling continuously?

%disable all controls except pause
set(handles.FirstFrame,'enable','off');
set(handles.PrevFrame,'enable','off');
set(handles.Play,'enable','off');
set(handles.Pause,'enable','on');
set(handles.NextFrame,'enable','off');
set(handles.LastFrame,'enable','off');

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
set(handles.ApertureFraction,'Enable','on');
set(handles.FrameRate,'Enable','on');

CurFrame = str2double(get(handles.CurFrame,'String'));
TotFrames = str2double(get(handles.TotFrames,'String'));

if (CurFrame > 1)
    set([handles.PrevFrame handles.FirstFrame],'Enable','on');
end

if (CurFrame < TotFrames)
    set([handles.NextFrame handles.LastFrame],'Enable','on');
end

set(handles.Pause,'Enable','off');

% --- Executes on button press in NextFrame.
function NextFrame_Callback(hObject, eventdata, handles)
% hObject    handle to NextFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%advance Frame
CurFrame = str2double(get(handles.CurFrame,'String')) + 1;
TotFrames = str2double(get(handles.NumFrames,'String'));

if (CurFrame > 1)
    set([handles.FirstFrame handles.PrevFrame],'enable','on');
else
    set([handles.FirstFrame handles.PrevFrame],'enable','off');
end

if (CurFrame == TotFrames)    
    set([handles.NextFrame handles.LastFrame],'enable','off');
else
    set([handles.NextFrame handles.LastFrame],'enable','on');
end

UpdateAnimationFrame(CurFrame, handles);

% --- Executes on button press in LastFrame.
function LastFrame_Callback(hObject, eventdata, handles)
% hObject    handle to LastFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%advance Frame
set([handles.FirstFrame handles.PrevFrame],'enable','on');
set([handles.NextFrame handles.LastFrame],'enable','off');

CurFrame = str2double(get(handles.NumFrames,'String'));
UpdateAnimationFrame(CurFrame, handles);


% --- Executes during object creation, after setting all properties.
function FastTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FastTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function NumFrames_Callback(hObject, eventdata, handles)
% hObject    handle to NumFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumFrames as text
%        str2double(get(hObject,'String')) returns contents of NumFrames as a double



function ApertureFraction_Callback(hObject, eventdata, handles)
% hObject    handle to ApertureFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ApertureFraction as text
%        str2double(get(hObject,'String')) returns contents of ApertureFraction as a double


% --- Executes during object creation, after setting all properties.
function ApertureFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ApertureFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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


% --- Executes on button press in CycleCheck.
function CycleCheck_Callback(hObject, eventdata, handles)
% hObject    handle to CycleCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CycleCheck


function UpdateAnimationFrame(CurFrame, handles)

if isnan(CurFrame)
    CurFrame = 1;
end
set(handles.CurFrame,'String',CurFrame);

NumFrames = str2double(get(handles.NumFrames,'String'));
ApertureFraction = str2double(get(handles.ApertureFraction,'String'));
if (get(handles.SlowTime,'Value'))
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


% --- Executes on button press in UnitsCheck.
function UnitsCheck_Callback(hObject, eventdata, handles)
% hObject    handle to UnitsCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UnitsCheck

%recompute sample spacing for full aperture
is_english = get(hObject,'Value');
set_res_fields(handles.meta.Grid.Col.SS, handles.CrossRangeSS, handles.CRSSUnits, is_english);
set_res_fields(handles.meta.Grid.Row.SS, handles.RangeSS, handles.SSUnits, is_english);

%change text on units combo boxes absed on mode and possibly units
if get(handles.Resolution,'Value')
    if get(hObject,'Value')
        %set as inches and ft
        set(handles.MinResUnits,'String','inches|feet');
        set(handles.MaxResUnits,'String','inches|feet');
        set(handles.StepSizeUnits,'String','inches|feet');
    else
        %set as cm and meters
        set(handles.MinResUnits,'String','cm|meters');
        set(handles.MaxResUnits,'String','cm|meters');
        set(handles.StepSizeUnits,'String','cm|meters');
    end
else
    %set as percent
    set(handles.MinResUnits,'String','Percent');
    set(handles.MaxResUnits,'String','Percent');
    set(handles.StepSizeUnits,'String','Percent');
end

SetResAnimationSettings(handles)

%this will update resolution
newpos(getPosition(handles.H), hObject);


% --- Executes during object creation, after setting all properties.
function MinMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function MinRes_Callback(hObject, eventdata, handles)
% hObject    handle to MinRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinRes as text
%        str2double(get(hObject,'String')) returns contents of MinRes as a double


% --- Executes on selection change in MinResUnits.
function MinResUnits_Callback(hObject, eventdata, handles)
% hObject    handle to MinResUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MinResUnits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MinResUnits


% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


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


% --- Executes on selection change in MaxResUnits.
function MaxResUnits_Callback(hObject, eventdata, handles)
% hObject    handle to MaxResUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MaxResUnits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MaxResUnits


% --- Executes during object creation, after setting all properties.
function text3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function StepSize_Callback(hObject, eventdata, handles)
% hObject    handle to StepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StepSize as text
%        str2double(get(hObject,'String')) returns contents of StepSize as a double


% --- Executes during object creation, after setting all properties.
function StepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in StepSizeUnits.
function StepSizeUnits_Callback(hObject, eventdata, handles)
% hObject    handle to StepSizeUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns StepSizeUnits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from StepSizeUnits


% --- Executes on button press in FullRangeBandwidth.
function FullRangeBandwidth_Callback(hObject, eventdata, handles)
% hObject    handle to FullRangeBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.FullAzBandwidth,'Enable','off');
else
    set(handles.FullAzBandwidth,'Enable','on');
end

%only modify imrect position if we are in ResMode
if get(handles.ResModeCheck,'Value')
    pos = getPosition(handles.H);
    if get(hObject,'Value')
        %set position to full bandwidth
        pos(2) = handles.ZpLimsRn(1);
        pos(4) = handles.ZpWidthRn;
    else
        pos = ForceCenterAperture(pos);
    end
    setPosition(handles.H,pos);
end
    
    
% --- Executes during object creation, after setting all properties.
function FullRangeBandwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FullRangeBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function CurFrameRes_Callback(~, eventdata, handles)
% hObject    handle to CurFrameRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurFrameRes as text
%        str2double(get(hObject,'String')) returns contents of CurFrameRes as a double


function TotFramesRes_Callback(hObject, eventdata, handles)
% hObject    handle to TotFramesRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TotFramesRes as text
%        str2double(get(hObject,'String')) returns contents of TotFramesRes as a double


% --- Executes on button press in FirstFrameRes.
function FirstFrameRes_Callback(hObject, eventdata, handles)
% hObject    handle to FirstFrameRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.CurFrameRes,'String',1);

set([handles.FirstFrameRes, handles.PrevFrameRes],'enable','off');
set([handles.NextFrameRes, handles.LastFrameRes],'enable','on');

UpdateResolutionFrame(handles);


% --- Executes on button press in PrevFrameRes.
function PrevFrameRes_Callback(hObject, eventdata, handles)
% hObject    handle to PrevFrameRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%go back Frame
CurFrame = str2double(get(handles.CurFrameRes,'String'));
if isnan(CurFrame)
    CurFrame = 1;
else
    CurFrame = CurFrame-1;
end
set(handles.CurFrameRes,'String',CurFrame);

if (CurFrame > 1)
    set([handles.FirstFrameRes, handles.PrevFrameRes],'enable','on');
else
    set([handles.FirstFrameRes, handles.PrevFrameRes],'enable','off');
end
set([handles.NextFrameRes, handles.LastFrameRes],'enable','on');

UpdateResolutionFrame(handles);


% --- Executes on button press in PlayRes.
function PlayRes_Callback(hObject, eventdata, handles)
% hObject    handle to PlayRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%play "resolution" movie

%Compute number of frames & validate inputs
MinRes = str2double(get(handles.MinRes,'String'));
MaxRes = str2double(get(handles.MaxRes,'String'));
StepSize = str2double(get(handles.StepSize,'String'));
MinResStrings = get(handles.MinResUnits,'String');
MaxResStrings = get(handles.MaxResUnits,'String');
StepSizeStrings = get(handles.StepSizeUnits,'String');
MinResValue = get(handles.MinResUnits,'Value');
MaxResValue = get(handles.MaxResUnits,'Value');
StepSizeValue = get(handles.StepSizeUnits,'Value');
MinResUnits = strtrim(MinResStrings(MinResValue,:));
MaxResUnits = strtrim(MaxResStrings(MaxResValue,:));
StepSizeUnits = strtrim(StepSizeStrings(StepSizeValue,:));

%convert to ft,meters or percent
if strcmp(MinResUnits,'inches')
    MinRes = MinRes/12;
elseif strcmp(MinResUnits,'cm')
    MinRes = MinRes/100;
end

if strcmp(MaxResUnits,'inches')
    MaxRes = MaxRes/12;
elseif strcmp(MaxResUnits,'cm')
    MaxRes = MaxRes/100;
end

if strcmp(StepSizeUnits,'inches')
    StepSize = StepSize/12;
elseif strcmp(StepSizeUnits,'cm')
    StepSize = StepSize/100;
end

if get(handles.Resolution,'Value')
    if (MaxRes>MinRes)
        msgbox('Max Resolution must be smaller than Min Resolution!!');
        return;
    end
    NumSteps = round((MinRes-MaxRes)/StepSize)+1;
else
    if (MaxRes<MinRes)
        msgbox('Max Resolution must be greater than Min Resolution!!');
        return;
    end
    NumSteps = round((MaxRes-MinRes)/StepSize)+1;
end

set(handles.CurFrameRes,'String',1);
set(handles.TotFramesRes,'String',NumSteps);

FrameRate = str2double(get(handles.FrameRate,'String'));

%disable all controls except pause
set([handles.FirstFrameRes,...
    handles.PrevFrameRes,...
    handles.PlayRes,...
    handles.NextFrameRes,...
    handles.LastFrameRes],'enable','off');
set(handles.PauseRes,'enable','on');

%determine if we are cycling continuously
CycleCheck = get(handles.CycleCheck,'Value');

handles.InLoop = 1;
guidata(hObject, handles);
CurFrame = 1;

while (handles.InLoop)
    %see if pause was pressed
    if (strcmp(get(handles.PlayRes,'Enable'),'on'))
        handles.InLoop = 0;
        break;
    end
    
    set(handles.CurFrameRes,'String',CurFrame);
    tic;
    
    %update display
    UpdateResolutionFrame(handles)
    
    %wait until time for next frame
    elapsedtime = toc;
    pausetime = 1/FrameRate - elapsedtime;                
    if (pausetime > 0)        
        pause(pausetime);
    else        
        pause(.05); %pause a little so the animation can be stopped
    end
    
    if (CurFrame == NumSteps)
        if (CycleCheck)
            CurFrame = 1;
        else
            %stop movie
            set(handles.PlayRes,'enable','on');
            set(handles.PauseRes,'enable','off');
            return;
        end
    else
        CurFrame = CurFrame + 1;
    end
end


% --- Executes on button press in PauseRes.
function PauseRes_Callback(hObject, eventdata, handles)
% hObject    handle to PauseRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.PlayRes,'Enable','on');

CurFrame = str2double(get(handles.CurFrameRes,'String'));
TotFrames = str2double(get(handles.TotFramesRes,'String'));

if (CurFrame > 1)
    set([handles.PrevFrameRes, handles.FirstFrameRes],'Enable','on');
end

if (CurFrame < TotFrames)
    set([handles.NextFrameRes, handles.LastFrameRes],'Enable','on');
end

set(handles.PauseRes,'Enable','off');


% --- Executes on button press in NextFrameRes.
function NextFrameRes_Callback(hObject, eventdata, handles)
% hObject    handle to NextFrameRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%advance Frame
CurFrame = str2double(get(handles.CurFrameRes,'String'));
TotFrames = str2double(get(handles.TotFramesRes,'String'));

if isnan(CurFrame)
    CurFrame = 1;
else
    CurFrame = CurFrame+1;
end

set(handles.CurFrameRes,'String',CurFrame);

if (CurFrame > 1)
    set([handles.FirstFrameRes, handles.PrevFrameRes],'enable','on');
else
    set([handles.FirstFrameRes, handles.PrevFrameRes],'enable','off');
end

if (CurFrame == TotFrames)
    set([handles.NextFrameRes, handles.LastFrameRes],'enable','off');
else
    set([handles.NextFrameRes, handles.LastFrameRes],'enable','on');
end

UpdateResolutionFrame(handles);


% --- Executes on button press in LastFrameRes.
function LastFrameRes_Callback(hObject, eventdata, handles)
% hObject    handle to LastFrameRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TotFrames = str2double(get(handles.TotFramesRes,'String'));
set(handles.CurFrameRes,'String',TotFrames);

set([handles.FirstFrameRes, handles.PrevFrameRes],'enable','on');
set([handles.NextFrameRes, handles.LastFrameRes],'enable','off');

UpdateResolutionFrame(handles);


% --- Executes during object creation, after setting all properties.
function MaxResUnits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxResUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ResModeCheck.
function ResModeCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ResModeCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ResModeCheck

if get(hObject,'Value')
    handles.ResMode = 1;
    %Enable/Disable Conrols
    set(handles.MaxRes,'enable','on');
    set(handles.MinRes,'enable','on');
    set(handles.StepSize,'enable','on');
    set(handles.PlayRes,'enable','on');
    set(handles.LastFrameRes,'enable','on');
    set(handles.FirstFrameRes,'enable','on');
    set(handles.NextFrameRes,'enable','off');
    set(handles.PrevFrameRes,'enable','off');
    setPositionConstraintFcn(handles.H,@ForceCenterAperture); 

    % Update handles structure
    guidata(hObject, handles);
else
    handles.ResMode = 0;
    %Enable/Disable Conrols
    set([handles.MaxRes,...
        handles.MinRes,...
        handles.StepSize,...
        handles.PlayRes,...
        handles.NextFrameRes,...
        handles.FirstFrameRes,...
        handles.PrevFrameRes,...
        handles.LastFrameRes],'enable','off');       
end

% Update handles structure
guidata(hObject, handles);


function pos = ForceCenterAperture(newpos)

%force the aperture to be centered and contain the same percentage in az
%and range.
handles = guidata(gca);

if get(handles.ResModeCheck,'Value')
    XPercent = newpos(3)/handles.ZpWidthAz;
    YPercent = newpos(4)/handles.ZpWidthRn;    
    
    if get(handles.GroundProjCheck,'Value') && strcmp(handles.meta.Grid.ImagePlane,'SLANT')
        %need to match resolutions, so we will match to the minimum ground
        %projected resolution (that way we don't have to check to see if we
        %have the resolution)
        AzRes = handles.meta.Grid.Col.ImpRespWid/cosd(handles.meta.SCPCOA.TwistAng)/XPercent;
        RnRes = handles.meta.Grid.Row.ImpRespWid/cosd(handles.meta.SCPCOA.GrazeAng)/YPercent;
        MinRes = max([AzRes RnRes]);
        
        AzPercent = (handles.meta.Grid.Col.ImpRespWid/cosd(handles.meta.SCPCOA.TwistAng))/MinRes;
        RnPercent = (handles.meta.Grid.Row.ImpRespWid/cosd(handles.meta.SCPCOA.GrazeAng))/MinRes;
    else
        %we will assume slant res is close to square. For percent aperture
        %we will just make them the same.
        %use minimum percent
        AzPercent = min([XPercent YPercent]);
        RnPercent = AzPercent;
    end
    %set new position
    pos(1) = round(handles.ZpWidthAz/2 - (AzPercent*handles.ZpWidthAz)/2) + ...
    handles.ZpLimsAz(1);
    pos(2) = round(handles.ZpWidthRn/2 - (RnPercent*handles.ZpWidthRn)/2) + ...
        handles.ZpLimsRn(1);
    pos(3) = round(handles.ZpWidthAz*AzPercent);
    pos(4) = round(handles.ZpWidthRn*RnPercent);
    
    %if full range/az BW is specified then set range to full
    if get(handles.FullRangeBandwidth,'Value')
        pos(2) = handles.ZpLimsRn(1);
        pos(4) = handles.ZpWidthRn;
    end
    if get(handles.FullAzBandwidth,'Value')
        pos(1) = handles.ZpLimsAz(1);
        pos(3) = handles.ZpWidthAz;
    end
else
    if verLessThan('matlab', '7.6')
        disp('Warning. Matlab version does not support imrect constraint function');
    else
        setPositionConstraintFcn(handles.H,...
            makeConstrainToRectFcn('imrect',handles.ZpLimsAz,handles.ZpLimsRn));
    end
    pos = newpos;
end


% --- Executes on button press in GroundProjCheck.
function GroundProjCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GroundProjCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SetResAnimationSettings(handles)


function SetResAnimationSettings(handles)

if get(handles.Resolution,'Value')
    % Computer all values in meters
    CrossRangeRes = handles.meta.Grid.Col.ImpRespWid;
    RangeRes = handles.meta.Grid.Row.ImpRespWid;        
    if get(handles.GroundProjCheck,'Value')
        CrossRangeRes = CrossRangeRes/cosd(handles.meta.SCPCOA.TwistAng);
        RangeRes = RangeRes/cosd(handles.meta.SCPCOA.GrazeAng);
    end
    ResMax = max([CrossRangeRes RangeRes]);
    ResMin = ResMax*10;
    StepSize = (ResMin-ResMax)/9;
    % Then display is request units.
    is_english = get(handles.UnitsCheck,'Value');
    set_res_fields(ResMax, handles.MaxRes, handles.MaxResUnits, is_english);
    set_res_fields(ResMin, handles.MinRes, handles.MinResUnits, is_english);
    set_res_fields(StepSize, handles.StepSize, handles.StepSizeUnits, is_english);
else
    %percent
    set(handles.MaxRes,'String',100);
    set(handles.MaxResUnits,'Value',1);
    set(handles.MinRes,'String',10);
    set(handles.MinResUnits,'Value',1);        
    set(handles.StepSize,'String',10);
    set(handles.StepSizeUnits,'Value',1);        
end


function ResModeChange(hObject, eventdata)

handles = guidata(gcf);

%change text on units combo boxes absed on mode and possibly units
Mode = get(handles.Resolution,'Value');
Units = get(handles.UnitsCheck,'Value');

if Mode
    if Units
        %set as inches and ft
        set(handles.MinResUnits,'String','inches|feet');
        set(handles.MaxResUnits,'String','inches|feet');
        set(handles.StepSizeUnits,'String','inches|feet');
    else
        %set as cm and meters
        set(handles.MinResUnits,'String','cm|meters');
        set(handles.MaxResUnits,'String','cm|meters');
        set(handles.StepSizeUnits,'String','cm|meters');
    end
    set(handles.GroundProjCheck,'enable','on');
else
    %set as percent
    set(handles.MinResUnits,'String','Percent');
    set(handles.MaxResUnits,'String','Percent');
    set(handles.StepSizeUnits,'String','Percent');
    set(handles.GroundProjCheck,'enable','off');
end

SetResAnimationSettings(handles)


function DomainChange(hObject, eventdata)

handles = guidata(gcf);

domain = get(handles.ImageRadio,'Value');

if domain
    %dispersed, only make JPEG save available
    set(handles.SICD,'Value',1);
    set(handles.SICD,'enable','on');
else
    %dispersed, only make JPEG save available
    set(handles.JPEG,'Value',1);
    set(handles.SICD,'enable','off');
end

function FilterChange(hObject, eventdata)

handles = guidata(gcf);
newpos(getPosition(handles.H), hObject);

function UpdateResolutionFrame(handles)

%Compute number of frames & validate inputs
MinRes = str2double(get(handles.MinRes,'String'));
MaxRes = str2double(get(handles.MaxRes,'String'));
StepSize = str2double(get(handles.StepSize,'String'));
MinResStrings = get(handles.MinResUnits,'String');
MaxResStrings = get(handles.MaxResUnits,'String');
StepSizeStrings = get(handles.StepSizeUnits,'String');
MinResValue = get(handles.MinResUnits,'Value');
MaxResValue = get(handles.MaxResUnits,'Value');
StepSizeValue = get(handles.StepSizeUnits,'Value');
MinResUnits = strtrim(MinResStrings(MinResValue,:));
MaxResUnits = strtrim(MaxResStrings(MaxResValue,:));
StepSizeUnits = strtrim(StepSizeStrings(StepSizeValue,:));

%convert to meters
if strcmp(MinResUnits,'inches')
    MinRes = (MinRes/12)*FEET_TO_METERS;
elseif strcmp(MinResUnits,'feet')
    MinRes = MinRes*FEET_TO_METERS;
elseif strcmp(MinResUnits,'cm')
    MinRes = MinRes/100;
end

if strcmp(MaxResUnits,'inches')
    MaxRes = (MaxRes/12)*FEET_TO_METERS;
elseif strcmp(MaxResUnits,'feet')
    MaxRes = MaxRes*FEET_TO_METERS;
elseif strcmp(MaxResUnits,'cm')
    MaxRes = MaxRes/100;
end

if strcmp(StepSizeUnits,'inches')
    StepSize = (StepSize/12)*FEET_TO_METERS;
elseif strcmp(StepSizeUnits,'feet')
    StepSize = StepSize*FEET_TO_METERS;
elseif strcmp(StepSizeUnits,'cm')
    StepSize = StepSize/100;
end

CurFrame = str2double(get(handles.CurFrameRes,'String'));

%determine percent or resolution for current step
MaxMin = get(handles.MaxMin,'Value');
if get(handles.Resolution,'Value')
    %resolution
    if MaxMin
        Res = MaxRes+(CurFrame-1)*StepSize;
    else
        Res = MinRes-(CurFrame-1)*StepSize;
    end
    
    %now determine necessary ap percent for resolution
    CrossRangeRes = handles.meta.Grid.Col.ImpRespWid;
    RangeRes = handles.meta.Grid.Row.ImpRespWid;

    if get(handles.GroundProjCheck,'Value')
        CrossRangeRes = CrossRangeRes/cosd(handles.meta.SCPCOA.TwistAng);
        RangeRes = RangeRes/cosd(handles.meta.SCPCOA.GrazeAng);
    end
    
    XPercent = CrossRangeRes/Res;
    YPercent = RangeRes/Res;
else
    %percent
    if MaxMin
        XPercent = (MaxRes-(CurFrame-1)*StepSize)/100;
    else
        XPercent = (MinRes+(CurFrame-1)*StepSize)/100;
    end
    YPercent = XPercent;
end

%set imrect position
pos(1) = round(handles.ZpWidthAz/2 - (XPercent*handles.ZpWidthAz)/2) + ...
    handles.ZpLimsAz(1);
pos(2) = round(handles.ZpWidthRn/2 - (YPercent*handles.ZpWidthRn)/2) + ...
    handles.ZpLimsRn(1);
pos(3) = round(handles.ZpWidthAz*XPercent);
pos(4) = round(handles.ZpWidthRn*YPercent);

%if full range/az BW is specified then set range to full
if get(handles.FullRangeBandwidth,'Value')
    pos(2) = handles.ZpLimsRn(1);
    pos(4) = handles.ZpWidthRn;
end
if get(handles.FullAzBandwidth,'Value')
    pos(1) = handles.ZpLimsAz(1);
    pos(3) = handles.ZpWidthAz;
end

%update imrect
setPosition(handles.H,pos);


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

if get(handles.JPEG,'Value')
    [fname, path] = uiputfile( {'*.jpg','jpg Files (*.jpg)'},'Save JPG File',pathstr);
    filename = strcat(path,fname);
    if get(handles.FullResCheck,'Value')
        %full res image of chip/phd
        if get(handles.ImageRadio,'Value')
            chip = abs(getFullResImage(handles));
        else
            chip = single(abs(handles.phasehistory));
            chipmean = mean(chip(:));
            chipstd = std(chip(:));
            chip(chip>(chipmean+3*chipstd)) = chipmean+3*chipstd;
            chip = chip./max(chip(:));
        end
        imwrite(chip,filename,'jpg');
    else
        %just write what is on the screen
        if get(handles.ImageRadio,'Value')
            frame = getframe(handles.image);
        else
            frame = getframe(handles.phd);
        end
        imwrite(frame.cdata,filename,'jpg');
    end   
    
else
    [fname, path] = uiputfile( {'*.ntf','SICD Files (*.ntf)'},'Save SICD File',pathstr);
    filename = strcat(path,fname);    
    
    %get complex image chip
    chip = getFullResImage(handles);
    chip = handles.fft_sp(fftshift(handles.fft_im((chip))));
    meta = handles.meta;
    meta.ImageFormation.Processing.Type = 'Subaperture';
    %TODO: update res, frequency support etc.
    meta.ImageData.NumRows = size(chip,1);
    meta.ImageData.NumCols = size(chip,2);
    meta = add_sicd_corners(meta); 
    
    writer = SICDWriter(filename,meta);
    writer.write_chip(chip.',[1 1]);
    delete(writer);    
end

msgbox('Image write complete');

setpref('matlab_sar_toolbox','last_used_directory',path); %store path


% --- Executes on button press in SaveMetaIcon.
function SaveMetaIcon_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMetaIcon (see GCBO)
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

[fname, path] = uiputfile( {'*.jpg','jpg Files (*.jpg)'},'Save JPG File',pathstr);
filename = strcat(path,fname);

frame = getframe(handles.metaicon);
imwrite(frame.cdata,filename,'jpg');

setpref('matlab_sar_toolbox','last_used_directory',path); %store path


% --- Executes on button press in Record.
function Record_Callback(hObject, eventdata, handles)
% hObject    handle to Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Recording = 1;
set(handles.Record,'enable','off');
set(handles.StopRecord,'enable','on');
handles.aviparams = [];

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in StopRecord.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to StopRecord (see GCBO)
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

set(handles.Record,'Enable','off');
set(handles.StopRecord,'Enable','off');

numframes = length(handles.aviparams);
set(handles.TotRecFrame,'String',numframes);

%get settings
FullRes = get(handles.FullRes,'Value');
ModeAnnotate = get(handles.ModeAnnotate,'Value');
FrameAnnotate = get(handles.FrameAnnotate,'Value');
ResAnnotate = get(handles.ResAnnotate,'Value');
AvgResAnnotate = get(handles.AvgResAnnotate,'Value');


framerate = str2double(get(handles.MovieFrameRate,'String'));

if FullRes
    [chip, handles] = getFullResImage(handles);
    [ny,nx,nz] = size(chip);

    %if chip is small then we will make bigger if annotation is selected,
    %otherwise it will not be legible
    MinSize = 1000;
    if ModeAnnotate || FrameAnnotate || ResAnnotate || AvgResAnnotate
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
        setPosition(handles.H,handles.aviparams(i).Pos);
        chip = abs(getFullResImage(handles));
        
%         %remap chip to specification
%         remaps = get(handles.remap_combo,'String');
%         remap = remaps{get(handles.remap_combo,'Value')};
%         chip = feval(remap,chip);        
       
        if ModeAnnotate || FrameAnnotate || ResAnnotate || AvgResAnnotate
            
            [ny,nx] = size(chip);
            if nx<MinSize && ny<MinSize
                if ny < nx
                    new_ny = MinSize;
                    new_nx = round((new_ny/ny)*nx);
                else
                    new_nx = MinSize;
                    new_ny = round((new_nx/nx)*ny);
                end
                chip = imresize(chip,[new_ny new_nx]);
                [ny,nx] = size(chip);
            end
            
            %set up annotation
            TextCount = 0;
            if ModeAnnotate
                TextCount = TextCount + 1;
                Annotation{TextCount} = sprintf('Mode: %s',handles.aviparams(i).Mode);
            end
            if FrameAnnotate
                TextCount = TextCount + 1;
                Annotation{TextCount} = sprintf('Frame %d of %d',i,numframes);
            end
            
            AzUnits = handles.aviparams(i).AzResUnits;
            RangeUnits = handles.aviparams(i).RangeResUnits;
            AzRes = handles.aviparams(i).AzRes;
            RangeRes = handles.aviparams(i).RangeRes;
            if (strcmp(AzUnits,'feet') || strcmp(AzUnits,'inches'))
                %english
                if strcmp(AzUnits,'inches')
                    AvgAzRes = AzRes/12;
                else
                    AvgAzRes = AzRes;
                end
                if strcmp(RangeUnits,'inches')
                    AvgRangeRes = RangeRes/12;
                else
                    AvgRangeRes = RangeRes;
                end
                AvgRes = (AvgAzRes+AvgRangeRes)/2;
                if (AvgRes < 1)
                    AvgRes = AvgRes*12;
                    Units = 'inches';
                else
                    Units = 'feet';
                end
            else
                %metric
                if strcmp(AzUnits,'cm')
                    AvgAzRes = AzRes/100;
                end
                if strcmp(RangeUnits,'cm')
                    AvgRangeRes = RangeRes/100;
                end
                AvgRes = (AvgAzRes+AvgRangeRes)/2;
                if (AvgRes < 1)
                    AvgRes = AvgRes*100;
                    Units = 'cm';
                else
                    Units = 'meters';
                end
            end
            
            if AvgResAnnotate
                TextCount = TextCount + 1;
                Annotation{TextCount} = sprintf('Avg Res: %2.1f %s',AvgRes,Units);
            else
                if  ResAnnotate
                    TextCount = TextCount + 1;
                    Annotation{TextCount} = sprintf('Az Res: %2.1f %s',AzRes,handles.aviparams(i).AzResUnits);
                    TextCount = TextCount + 1;
                    Annotation{TextCount} = sprintf('Range Res: %2.1f %s',RangeRes,handles.aviparams(i).RangeResUnits);
                end
            end
        end
    else
        %Screen Resolution, just update imrect pos and let callback update
        setPosition(handles.H,handles.aviparams(i).Pos);
        %set CLim
        set(handles.image,'CLim',handles.aviparams(i).CLim);
        %set scroll panel to selected view
        handles.apiSP1.setMagnificationAndCenter(handles.aviparams(i).Mag,...
        handles.aviparams(i).Center(1),handles.aviparams(i).Center(2));
        axes(handles.image);
        F = getframe( handles.image );
        chip = F.cdata(:,:,:);
        if strcmpi(polstring,'amplitude')
            animated(:,:,1,i) = F.cdata(:,:,1);
        else
            if i == 1
                [animated, cmap] = rgb2ind(F.cdata, 256, 'nodither');
            else
                animated(:,:,1,i) = rgb2ind(F.cdata, cmap, 'nodither');
            end
        end
        TextCount = 0; % Show no text annotations with screen resolution
    end

    %add annotation to chip if specified
    if TextCount > 0
        ax = round(nx/4);
        ay = round(ny*TextCount*0.02);
        aimage = ones(ay,ax)*255;
        rowheight = floor(ay/TextCount);
        spacing = round(rowheight*0.1);
        FontSize = rowheight - 2*spacing;
        for j=1:TextCount
            startpos = spacing +(j-1)*rowheight;
            aimage = AddTextToImage(aimage,Annotation{j},[startpos spacing],0,'Arial',FontSize);
        end
        
        if strcmpi(polstring,'amplitude')
            chip(ny-ay+1:end,1:ax) = aimage;
        else
            chip(ny-ay+1:end,1:ax,1) = aimage;
            chip(ny-ay+1:end,1:ax,2) = aimage;
            chip(ny-ay+1:end,1:ax,3) = aimage;
        end
    end
    
    %store chip frame for animated GIF
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
    set(handles.RecFrame,'String',i);
end

DelayTime = 1/framerate;
if strcmpi(polstring,'amplitude')
    imwrite(animated, filename, 'DelayTime', DelayTime, 'LoopCount', inf);
else
    imwrite(animated, cmap, filename, 'DelayTime', DelayTime, 'LoopCount', inf);
end

set(handles.TotRecFrame,'String','');
set(handles.RecFrame,'String','');

set(handles.Record,'Enable','on');
set(handles.StopRecord,'Enable','off');

handles.Recording = 0;

% Update handles structure
guidata(hObject, handles);


function RecFrame_Callback(hObject, eventdata, handles)
% hObject    handle to RecFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RecFrame as text
%        str2double(get(hObject,'String')) returns contents of RecFrame as a double


% --- Executes during object creation, after setting all properties.
function RecFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RecFrame (see GCBO)
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


% --- Executes on button press in FrameAnnotate.
function FrameAnnotate_Callback(hObject, eventdata, handles)
% hObject    handle to FrameAnnotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FrameAnnotate


% --- Executes on button press in ResAnnotate.
function ResAnnotate_Callback(hObject, eventdata, handles)
% hObject    handle to ResAnnotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ResAnnotate


function [chip,handles] = getFullResImage(handles)

%return full res complex image chip for current position
pos = getPosition(handles.H);
xmin = max(handles.ZpLimsAz(1),floor(pos(1)));
xmax = min(handles.ZpLimsAz(2),floor(pos(1)+pos(3)-1));
ymin = max(handles.ZpLimsRn(1),floor(pos(2)));
ymax = min(handles.ZpLimsRn(2),floor(pos(2)+pos(4)-1));

%need to modify xposition if flight is left
[ny, nx, nz] = size(handles.phasehistory);
if isfield(handles.meta,'SCPCOA') && isfield(handles.meta.SCPCOA,'SideOfTrack') && ...
        strcmp(handles.meta.SCPCOA.SideOfTrack,'L')
    xminold = xmin;
    xmin = round(nx-xmax);
    if (xmin<1);xmin=1;end;
    xmax = round(nx-xminold);    
end

%get filter type
if get(handles.FilterNone,'Value')
    WindowFilter = 0;
elseif get(handles.FilterGaussian,'Value')
    WindowFilter = 1;
elseif get(handles.Filterx4,'Value')
    WindowFilter = 2;
elseif get(handles.FilterHamming,'Value')
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
    phdchip = FilterPHDChip(phdchip,WindowFilter,FilterRange,FilterAz); 
    Ap(ymin:ymax,xmin:xmax,1:size(phdchip,3)) = phdchip;    
end

[chip,handles] = makeDisplayable(handles,Ap);

%chip = handles.fft_im(Ap);


% --- Executes on button press in FullResCheck.
function FullResCheck_Callback(hObject, eventdata, handles)
% hObject    handle to FullResCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FullResCheck


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


% --- Executes on button press in AvgResAnnotate.
function AvgResAnnotate_Callback(hObject, eventdata, handles)
% hObject    handle to AvgResAnnotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AvgResAnnotate


% --- Executes on button press in ModeAnnotate.
function ModeAnnotate_Callback(hObject, eventdata, handles)
% hObject    handle to ModeAnnotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ModeAnnotate


function CrossRangeGroundRes_Callback(hObject, eventdata, handles)
% hObject    handle to CrossRangeGroundRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossRangeGroundRes as text
%        str2double(get(hObject,'String')) returns contents of CrossRangeGroundRes as a double


% --- Executes during object creation, after setting all properties.
function CrossRangeGroundRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossRangeGroundRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RangeGroundRes_Callback(hObject, eventdata, handles)
% hObject    handle to RangeGroundRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RangeGroundRes as text
%        str2double(get(hObject,'String')) returns contents of RangeGroundRes as a double


% --- Executes during object creation, after setting all properties.
function RangeGroundRes_CreateFcn(hObject, ~, handles)
% hObject    handle to RangeGroundRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
    
newpos(getPosition(handles.H), hObject);


% --- Executes on button press in FullAzBandwidth.
function FullAzBandwidth_Callback(hObject, eventdata, handles)
% hObject    handle to FullAzBandwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    set(handles.FullRangeBandwidth,'Enable','off');
else
    set(handles.FullRangeBandwidth,'Enable','on');
end

%only modify imrect position if we are in ResMode
if get(handles.ResModeCheck,'Value')
    pos = getPosition(handles.H);
    if get(hObject,'Value')
        %set position to full bandwidth
        pos(1) = handles.ZpLimsAz(1);
        pos(3) = handles.ZpWidthAz;
    else
        pos = ForceCenterAperture(pos);
    end
    setPosition(handles.H,pos);
end


% --- Executes during object creation, after setting all properties.
function FullRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FullRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not c


% --- Executes during object creation, after setting all properties.
function TotFramesRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotFramesRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function CurFrameRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurFrameRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function MinRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function StepSizeUnits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepSizeUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on button press in SlowTime.
function SlowTime_Callback(hObject, eventdata, handles)
% hObject    handle to SlowTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SlowTime


% --- Executes on button press in FastTime.
function FastTime_Callback(hObject, eventdata, handles)
% hObject    handle to FastTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FastTime


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


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

set(parent.CenterPHDCheckAz,'Value',false);
set(parent.CenterPHDCheckRng,'Value',false);
process_phd(handles.parent,parent);
close(get(hObject,'Parent'));


% --- Executes on selection change in remap_combo.
function remap_combo_Callback(hObject, eventdata, handles)
% hObject    handle to remap_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')}, 'linearremap')
    set(handles.ImageContrast,'Visible','on');
else
    set(handles.ImageContrast,'Visible','off');
end
newpos(getPosition(handles.H),hObject); % Apply new remap


% --- Executes during object creation, after setting all properties.
function remap_combo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to remap_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in View3DVolume.
function View3DVolume_Callback(hObject, eventdata, handles)
% hObject    handle to View3DVolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get animation settings 
NumFrames = str2double(get(handles.NumFrames,'String'));
ApFrac = str2double(get(handles.ApertureFraction,'String'));

SlowFlag = get(handles.SlowTime,'Value');

%get filter type
if get(handles.FilterNone,'Value')
    WindowFilter = 0;
elseif get(handles.FilterGaussian,'Value')
    WindowFilter = 1;
elseif get(handles.Filterx4,'Value')
    WindowFilter = 2;
elseif get(handles.FilterHamming,'Value')
    WindowFilter = 3;
else
    WindowFilter = 4;
end

CollectTime = handles.meta.Timeline.CollectDuration;

%get display settings
settings.AlphaFactor = str2double(get(handles.AlphaFactor,'String'));
settings.Percentile = str2double(get(handles.Percentile,'String'));
settings.CMapCutoff = str2double(get(handles.CMapCutoff,'String'));

%display 3D Volume
if SlowFlag
    Bounds = [0 handles.meta.Timeline.CollectDuration];
else
    CenterFreq = handles.meta.Grid.Row.KCtr * SPEED_OF_LIGHT( ) / 2;
    BW = handles.meta.Grid.Row.ImpRespBW * SPEED_OF_LIGHT( ) / 2;
    Bounds = [(CenterFreq-BW/2)/1e6 (CenterFreq+BW/2)/1e6];
end

DisplaySubApVolume(handles.phasehistory,handles.ZpLimsAz,handles.ZpLimsRn,NumFrames,ApFrac,SlowFlag,WindowFilter,Bounds,settings);


function AlphaFactor_Callback(hObject, eventdata, handles)
% hObject    handle to AlphaFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AlphaFactor as text
%        str2double(get(hObject,'String')) returns contents of AlphaFactor as a double


% --- Executes during object creation, after setting all properties.
function AlphaFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CMapCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to CMapCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CMapCutoff as text
%        str2double(get(hObject,'String')) returns contents of CMapCutoff as a double


% --- Executes during object creation, after setting all properties.
function CMapCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CMapCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Percentile_Callback(hObject, eventdata, handles)
% hObject    handle to Percentile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Percentile as text
%        str2double(get(hObject,'String')) returns contents of Percentile as a double


% --- Executes during object creation, after setting all properties.
function Percentile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Percentile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AlwaysApply.
function AlwaysApply_Callback(hObject, eventdata, handles)
% hObject    handle to AlwaysApply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newpos(getPosition(handles.H), hObject);


% --- Executes on selection change in image_rep_combo.
function image_rep_combo_Callback(hObject, eventdata, handles)
% hObject    handle to image_rep_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_rep_combo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_rep_combo

newpos(getPosition(handles.H), hObject);


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


% --- Executes on button press in polar_vv_togglebutton.
function polar_vv_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to polar_vv_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polar_vv_togglebutton

if isfield(handles,'TxRcvPolarizationProc')
    handles.polar_PHD_idx = find(strcmpi('V:V',handles.TxRcvPolarizationProc));
end

updatePolarToggle(handles)
process_phd(hObject,handles)
newpos(getPosition(handles.H), hObject);



% --- Executes on button press in polar_vh_togglebutton.
function polar_vh_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to polar_vh_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polar_vh_togglebutton

if isfield(handles,'TxRcvPolarizationProc')
    handles.polar_PHD_idx = find(strcmpi('V:H',handles.TxRcvPolarizationProc));
end

updatePolarToggle(handles)
process_phd(hObject,handles)
newpos(getPosition(handles.H), hObject);


% --- Executes on button press in polar_hh_togglebutton.
function polar_hh_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to polar_hh_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polar_hh_togglebutton

if isfield(handles,'TxRcvPolarizationProc')
    handles.polar_PHD_idx = find(strcmpi('H:H',handles.TxRcvPolarizationProc));
end

updatePolarToggle(handles)
process_phd(hObject,handles)
newpos(getPosition(handles.H), hObject);


% --- Executes on button press in polar_hv_togglebutton.
function polar_hv_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to polar_hv_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of polar_hv_togglebutton

if isfield(handles,'TxRcvPolarizationProc')
    handles.polar_PHD_idx = find(strcmpi('H:V',handles.TxRcvPolarizationProc));
end

updatePolarToggle(handles)
process_phd(hObject,handles)
newpos(getPosition(handles.H), hObject);



function updatePolarToggle(handles)
% updatePolarToggle updated toggles for polarity
%
% updatePolarToggle(handles)
%
% INPUTS:
%   handles: figure handles
%
% Created by Pat Cutler January 2014 (NGA)

% reset all polor toggles
set(handles.polar_vv_togglebutton,'enable','on','value',0)
set(handles.polar_vh_togglebutton,'enable','on','value',0)
set(handles.polar_hh_togglebutton,'enable','on','value',0)
set(handles.polar_hv_togglebutton,'enable','on','value',0)

if isfield(handles,'TxRcvPolarizationProc')
    %disable toggles for unavailable
    if ~any(strcmpi('V:V',handles.TxRcvPolarizationProc))
        set(handles.polar_vv_togglebutton,'enable','off')
    end
    if ~any(strcmpi('V:H',handles.TxRcvPolarizationProc))
        set(handles.polar_vh_togglebutton,'enable','off')
    end
    if ~any(strcmpi('H:H',handles.TxRcvPolarizationProc))
        set(handles.polar_hh_togglebutton,'enable','off')
    end
    if ~any(strcmpi('H:V',handles.TxRcvPolarizationProc))
        set(handles.polar_hv_togglebutton,'enable','off')
    end
    
    %change value of toggle for PHD polarization
    switch handles.TxRcvPolarizationProc{handles.polar_PHD_idx}
        case 'V:V'
            set(handles.polar_vv_togglebutton,'value',1)
        case 'V:H'
            set(handles.polar_vh_togglebutton,'value',1)
        case 'H:H'
            set(handles.polar_hh_togglebutton,'value',1)
        case 'H:V'
            set(handles.polar_hv_togglebutton,'value',1)
    end
    
else
    %no polarization information available
    set(handles.polar_vv_togglebutton,'enable','off','value',0)
    set(handles.polar_vh_togglebutton,'enable','off','value',0)
    set(handles.polar_hh_togglebutton,'enable','off','value',0)
    set(handles.polar_hv_togglebutton,'enable','off','value',0)
end

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
    api_handle.addNewPositionCallback(@newlinepos);
    set(handles.Measure,'String','Remove');
    handles = ComputeDistance(handles);
else
    delete(handles.MeasureLine);
    handles.MeasureLine = [];
    set(handles.Measure,'String','Measure');
    set(handles.Distance,'String','');
end

guidata(hObject, handles);

function newlinepos(line)

handles = guidata(gcbo);
ComputeDistance(handles);

function handles = ComputeDistance(handles)

meta = handles.meta;

%compute line distance and heading in the ground plane
pos = getPosition(handles.MeasureLine);
StartPos = pos(1,:);
StopPos = pos(2,:);

DeltaX = (StartPos(1)-StopPos(1))*meta.Grid.Col.SS/cosd(meta.SCPCOA.TwistAng);
DeltaY = (StartPos(2)-StopPos(2))*meta.Grid.Row.SS/cosd(meta.SCPCOA.GrazeAng);

Distance = sqrt(DeltaX*DeltaX+DeltaY*DeltaY);

%alternate method...get lat/lon for each point and then compute the
%distance between the points
lla1 = point_slant_to_ground(StartPos',meta);
lla2 = point_slant_to_ground(StopPos',meta);

[d1km,d2km]=lldistkm(lla1(1:2),lla2(1:2));

Distance = d1km*1000;

set(handles.Distance,'String',sprintf('%.2f',Distance));

function Distance_Callback(hObject, eventdata, handles)
% hObject    handle to Distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Distance as text
%        str2double(get(hObject,'String')) returns contents of Distance as a double


% --- Executes during object creation, after setting all properties.
function Distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Distance (see GCBO)
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