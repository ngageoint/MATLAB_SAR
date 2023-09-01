function varargout = PulseExplorer(varargin)
% PULSEEXPLORER MATLAB code for PulseExplorer.fig
%      PULSEEXPLORER, by itself, creates a new PULSEEXPLORER or raises the existing
%      singleton*.
%
%      H = PULSEEXPLORER returns the handle to a new PULSEEXPLORER or the handle to
%      the existing singleton*.
%
%      PULSEEXPLORER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PULSEEXPLORER.M with the given input arguments.
%
%      PULSEEXPLORER('Property','Value',...) creates a new PULSEEXPLORER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PulseExplorer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PulseExplorer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PulseExplorer

% Last Modified by GUIDE v2.5 31-Jan-2014 13:27:28

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PulseExplorer_OpeningFcn, ...
                   'gui_OutputFcn',  @PulseExplorer_OutputFcn, ...
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


% --- Executes just before PulseExplorer is made visible.
function PulseExplorer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PulseExplorer (see VARARGIN)

% Choose default command line output for PulseExplorer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PulseExplorer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%set defaults
set(handles.FrameRate,'String',5);
set(handles.PulseIncrement,'String',1);

handles.Recording = 0;

%initialize metaicon holder
set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');
set(handles.metaicon,'Color',[51/255 102/255 153/255]);

%disable appropriate controls (off until data is actually loaded)
set(handles.CurPulse,'Enable','off');
set(handles.TotPulses,'Enable','off');
set(handles.FirstPulse,'Enable','off');
set(handles.PrevPulse,'Enable','off');
set(handles.PausePulse,'Enable','off');
set(handles.PlayPulse,'Enable','off');
set(handles.NextPulse,'Enable','off');
set(handles.LastPulse,'Enable','off');
set(handles.StartRecording,'Enable','off');
set(handles.SnapImage,'Enable','off');

%load pulse display names (these are m files in the PulseDisplay folder)
path=fileparts(mfilename('fullpath'));
fullfilelist=dir(fullfile(path, 'PulseDisplay', 'pe_disp_*.m'));
filelist=regexprep({fullfilelist.name}, 'pe_disp_(.*)\.m', '$1');
set(handles.PulseDisplayCombo,'String',filelist);

%set combo box to Reramped if available
foo = strcmpi('RFSignal',filelist);
if max(foo) > 0
    [val index] = max(foo);
else
    %just set to first item
    index = 1;
end
set(handles.PulseDisplayCombo,'Value',index);

%handles input filename if specified
p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.parse(varargin{:});
if length(p.Results.filename) > 0
    set(handles.filename,'String',p.Results.filename);
    handles = LoadFile(p.Results.filename,handles);    
end

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = PulseExplorer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Contrast.
function Contrast_Callback(hObject, eventdata, handles)
% hObject    handle to Contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


axes_handles = findobj(handles.display_panel,'Type','Axes');
if exist('imcontrast') %#ok<EXIST>
    h = imcontrast(axes_handles(1)); % Adjust a single axes
else
    h = imcontrast2(axes_handles(1)); % Adjust a single axes
end
    
% Set all others to the same contrast limits
set(h,'DeleteFcn',@(hObject, eventdata, handles) set(axes_handles,'CLim',get(axes_handles(1),'CLim')));

% --- Executes on selection change in ColormapCombo.
function ColormapCombo_Callback(hObject, eventdata, handles)
% hObject    handle to ColormapCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ColormapCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ColormapCombo

val = get(hObject,'Value');
strings = get(hObject,'String');
colormap(get(handles.figure1,'CurrentAxes'), strings{val});
% Each figure has only a single colormap, so setting the colormap for one
% axes sets them all.


% --- Executes during object creation, after setting all properties.
function ColormapCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColormapCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurPulse_Callback(hObject, eventdata, handles)
% hObject    handle to CurPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurPulse as text
%        str2double(get(hObject,'String')) returns contents of CurPulse as a double

handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function CurPulse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TotPulses_Callback(hObject, eventdata, handles)
% hObject    handle to TotPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TotPulses as text
%        str2double(get(hObject,'String')) returns contents of TotPulses as a double


% --- Executes during object creation, after setting all properties.
function TotPulses_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PulseDisplayCombo.
function PulseDisplayCombo_Callback(hObject, eventdata, handles)
% hObject    handle to PulseDisplayCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PulseDisplayCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PulseDisplayCombo

handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PulseDisplayCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PulseDisplayCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FirstPulse.
function FirstPulse_Callback(hObject, eventdata, handles)
% hObject    handle to FirstPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.CurPulse,'String',1);
handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PrevPulse.
function PrevPulse_Callback(hObject, eventdata, handles)
% hObject    handle to PrevPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurPulse = str2double(get(handles.CurPulse,'String'));
Increment = str2double(get(handles.PulseIncrement,'String'));
CurPulse = CurPulse-Increment;
set(handles.CurPulse,'String',CurPulse);
handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PausePulse.
function PausePulse_Callback(hObject, eventdata, handles)
% hObject    handle to PausePulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.PlayPulse,'Enable','on');
set(handles.PausePulse,'Enable','off');

% --- Executes on button press in PlayPulse.
function PlayPulse_Callback(hObject, eventdata, handles)
% hObject    handle to PlayPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.InLoop = 1;
set(handles.PausePulse,'Enable','on');
set(handles.PlayPulse,'Enable','off');

FrameRate = str2double(get(handles.FrameRate,'String'));
FrameSpace = 1/FrameRate;

while (handles.InLoop)
    
    %see if pause was pressed
    if (strcmp(get(handles.PlayPulse,'Enable'),'on'))
        handles.InLoop = 0;
        break;
    end
    
    localtic = tic;
    
    CurPulse = str2double(get(handles.CurPulse,'String'));
    Increment = str2double(get(handles.PulseIncrement,'String'));
    CurPulse = CurPulse+Increment;
    set(handles.CurPulse,'String',CurPulse);
    handles = UpdatePulseDisplay(handles);
    
    % Update handles structure
    guidata(hObject, handles);
    
    elapsed = toc(localtic);
    
    if (elapsed < FrameSpace)
        pause(FrameSpace-elapsed);
    else
        %pause a little so we can interrupt movie
        pause(.05);
    end
end


% --- Executes on button press in NextPulse.
function NextPulse_Callback(hObject, eventdata, handles)
% hObject    handle to NextPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurPulse = str2double(get(handles.CurPulse,'String'));
Increment = str2double(get(handles.PulseIncrement,'String'));
CurPulse = CurPulse + Increment;
set(handles.CurPulse,'String',CurPulse);
handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in LastPulse.
function LastPulse_Callback(hObject, eventdata, handles)
% hObject    handle to LastPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NumPulses = double(handles.meta.Data.Channel(1).NumVectors);
set(handles.CurPulse,'String',NumPulses);
handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function PulseSlider_Callback(hObject, eventdata, handles)
% hObject    handle to PulseSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

NumPulses = double(handles.meta.Data.Channel(1).NumVectors);
CurPulse = round(get(handles.PulseSlider,'Value')*NumPulses);
set(handles.CurPulse,'String',CurPulse);
handles = UpdatePulseDisplay(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PulseSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PulseSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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



function PulseIncrement_Callback(hObject, eventdata, handles)
% hObject    handle to PulseIncrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PulseIncrement as text
%        str2double(get(hObject,'String')) returns contents of PulseIncrement as a double


% --- Executes during object creation, after setting all properties.
function PulseIncrement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PulseIncrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StartRecording.
function StartRecording_Callback(hObject, eventdata, handles)
% hObject    handle to StartRecording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmpi(get(handles.StartRecording,'String'),'Start Recording')
    handles.Recording = 1;
    set(handles.StartRecording,'String','Stop Recording');
    handles.FrameList = [];
else
    %pause animation if playing
    set(handles.PlayPulse,'Enable','on');
    set(handles.PausePulse,'Enable','off');
    
    handles.Recording = 0;
    
    %get output filename (this will be an animated GIF)
    %load last path
    if ispref('matlab_sar_toolbox','last_used_directory')
        startpath = getpref('matlab_sar_toolbox','last_used_directory');
        if ~ischar(startpath)||~exist(startpath,'dir')
            startpath = pwd;
        end
    else
        startpath = pwd;
    end

    %get jpg file name 
    [fname,pathstr] = uiputfile('*.gif','Save Animated GIF',startpath);
    outfile = [pathstr fname];
     
    numframes = length(handles.FrameList);
    if numframes <=0
        return;
    end
    animated(1,1,1,numframes) = 0;
    FrameRate = str2double(get(handles.FrameRate,'String'));
    %loop through saved frames and capture each frame
    for i=1:length(handles.FrameList)
        set(handles.CurPulse,'String',handles.FrameList(i));
        handles = UpdatePulseDisplay(handles);
        
        %Trim Frame
        F = getframe(gcf);
        [ny nx nz] = size(F.cdata);
        StartX = round(0.01*nx)+1;
        StopX = nx - round(0.01*nx);
        StartY = round(0.067*ny)+1;
        StopY = round(0.79*ny);
      
        if (i==1)
            [animated, cmap] = rgb2ind(F.cdata(StartY:StopY,StartX:StopX,:), 256, 'nodither');
        else
            animated(:,:,1,i) = rgb2ind(F.cdata(StartY:StopY,StartX:StopX,:), cmap, 'nodither');
        end
    end
    
    DelayTime = 1/FrameRate;
    imwrite(animated, cmap, outfile, 'DelayTime', DelayTime, ...
    'LoopCount', inf);
    
    set(handles.StartRecording,'String','Start Recording');
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SnapImage.
function SnapImage_Callback(hObject, eventdata, handles)
% hObject    handle to SnapImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    startpath = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(startpath)||~exist(startpath,'dir')
        startpath = pwd;
    end
else
    startpath = pwd;
end

%get jpg file name 
[fname,pathstr] = uiputfile('*.jpg','Save Image File',startpath);
outfile = [pathstr fname];

%get frame
F = getframe(gcf);
[ny nx nz] = size(F.cdata);
StartX = round(0.01*nx)+1;
StopX = nx - round(0.01*nx);
StartY = round(0.067*ny)+1;
StopY = round(0.79*ny);

imwrite(F.cdata(StartY:StopY,StartX:StopX,:),outfile,'jpg');

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

%get phase history
[fname, pathstr] = uigetfile( sar_file_extensions('phd'),...
    'Open Phase History File', pathstr);

if (isnumeric(fname))
    return;
end
file = strcat(pathstr,fname);
set(handles.filename,'String',file);

handles = LoadFile(file,handles);

%store path
setpref('matlab_sar_toolbox','last_used_directory',pathstr);

% Update handles structure
guidata(hObject, handles);

function handles = LoadFile(file,handles)

delete(get(handles.display_panel,'Children'));

if isfield(handles,'ph_reader')
    handles.ph_reader.close();
end

handles.phname = file;

%open reader
handles.ph_reader = open_ph_reader(file);
meta = handles.ph_reader.get_meta();
handles.meta = meta;

%set Channel and Display Selector
path=fileparts(mfilename('fullpath'));
if isfield(handles.meta.Data,'NumCRSDChannels')  % Check raw first if available
    NumChannels = handles.meta.Data.NumCRSDChannels;
    fullfilelist=dir(fullfile(path, 'PulseDisplay', 'pe_disp_*.m'));
elseif isfield(handles.meta.Data,'NumCPHDChannels')
    if ~strcmp(meta.Global.DomainType, 'FX')
        error('PULSEEXPLORER:NOT_FX_DOMAIN','Unable to process TOA domain data.');
    end
    NumChannels = handles.meta.Data.NumCPHDChannels;
    % Set available display types
    fullfilelist=dir(fullfile(path, 'PulseDisplay', 'CPHD', 'pe_disp_*.m'));
end
% Set available display types
filelist=regexprep({fullfilelist.name}, 'pe_disp_(.*)\.m', '$1');
set(handles.PulseDisplayCombo,'String',filelist);
index = find(strcmpi('CPHD',filelist));
if isempty(index), index = 1; end
set(handles.PulseDisplayCombo,'Value',index);
% Channel
ChannelStrings=cell(1,NumChannels);
if isfield(meta,'Channel') && isfield(meta.Channel,'Parameters')
    for i=1:NumChannels
        if isfield(meta.Channel.Parameters(i),'Identifier')
            ChannelStrings{i} = meta.Channel.Parameters(i).Identifier;
        elseif isfield(meta.Channel.Parameters(i),'RcvPol') && ...
                isfield(meta.Channel.Parameters(i),'SARImaging') && ...
                isfield(meta.Channel.Parameters(i).SARImaging,'TxPol')
            ChannelStrings{i} =  sprintf('Channel %d: %s',i, ...
                [meta.Channel.Parameters(i).RcvPol ':' ...
                meta.Channel.Parameters(i).SARImaging.TxPol]);
        else
            ChannelStrings{i} =  sprintf('Channel: %d', i);
        end
    end
end

set(handles.ChannelCombo,'String',ChannelStrings);
set(handles.ChannelCombo,'Value',1);

%set CurPulse and NumPulses
set(handles.CurPulse,'String',500);
set(handles.PulseSlider,'Value',500/handles.meta.Data.Channel(1).NumVectors);
set(handles.TotPulses,'String',handles.meta.Data.Channel(1).NumVectors);

%draw metaicon
MetaIcon(file, 'handle', handles.metaicon);

%update display
handles = UpdatePulseDisplay(handles);

%enable appropriate controls
set(handles.CurPulse,'Enable','on');
set(handles.PlayPulse,'Enable','on');
set(handles.NextPulse,'Enable','on');
set(handles.LastPulse,'Enable','on');
set(handles.PrevPulse,'Enable','on');
set(handles.FirstPulse,'Enable','on');
set(handles.SnapImage,'Enable','on');
set(handles.StartRecording,'Enable','on');


function handles = UpdatePulseDisplay(handles)

%get function name to call
strings = get(handles.PulseDisplayCombo,'String');
val = get(handles.PulseDisplayCombo,'Value');
pulsedisplay = ['pe_disp_' strings{val}];

%get pulse number and channel
PulseNum = str2double(get(handles.CurPulse,'String'));
Channel = get(handles.ChannelCombo,'Value');

% Requested pulse number might be out of range.  If so, must fix.
if isfield(handles.meta,'extra')
    min_pulse = 1 - handles.meta.extra.nonimaging_pulses_pre;
else
    min_pulse = 1;
end

if (PulseNum <= min_pulse);
    PulseNum = min_pulse;
    set(handles.PrevPulse,'Enable','off');
else
    set(handles.PrevPulse,'Enable','on');
end

if (PulseNum >= double(handles.meta.Data.Channel(Channel).NumVectors))
    PulseNum = double(handles.meta.Data.Channel(Channel).NumVectors);
    set(handles.NextPulse, 'Enable','off');
    set(handles.LastPulse, 'Enable','off');
else
    set(handles.NextPulse, 'Enable','on');
    set(handles.LastPulse, 'Enable','on');
end

if (PulseNum == 1)
    set(handles.FirstPulse,'Enable','off');
else
    set(handles.FirstPulse,'Enable','on');
end

set(handles.CurPulse,'String',PulseNum); % Must adjust if out-of-range input was given
set(handles.PulseSlider,'Value',min(1,max(0,double(PulseNum)/double(handles.meta.Data.Channel(Channel).NumVectors))));

%store pulse number if we are recording a movie
if (handles.Recording == 1)
    handles.FrameList = horzcat(handles.FrameList,PulseNum);
end

%get CLim from existing plots
axes_handles = findobj(handles.display_panel, 'Type', 'axes');
if ~isempty(axes_handles)
    CLim = get(axes_handles(1),'CLim');  
else
    CLim = [0 0];
end

%now delete any old plots
delete(get(handles.display_panel,'Children'));

%call specified pulse display function
feval(pulsedisplay,handles.display_panel,handles.ph_reader,PulseNum,Channel,CLim,handles.phname);
datacursormode(handles.figure1, 'on'); % Allow user to query data by position


% --- Executes on selection change in ChannelCombo.
function ChannelCombo_Callback(hObject, eventdata, handles)
% hObject    handle to ChannelCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChannelCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChannelCombo

handles = UpdatePulseDisplay(handles);


% --- Executes during object creation, after setting all properties.
function ChannelCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Disp_Metadata.
function Disp_Metadata_Callback(hObject, eventdata, handles)
% hObject    handle to Disp_Metadata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'meta')
    MetaViewer(handles.meta);
else
    MetaViewer;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////