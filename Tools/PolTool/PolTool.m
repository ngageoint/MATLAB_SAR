function varargout = PolTool(varargin)
%PolTool Interactive Polarization Combination Tool
%
% Applies weighted complex combination of linear polarized data.  For
% Quad-Pol this allows us to synthesize any transmit and recieve
% polarimetric configuration.  Dual Pol allows us to only synthesize one 
% side (typically Rx). 
% Polarization is defined by angle and ellipticitycontrol.  Angle is defined from
% 0-360 (where 0 is H).  EllipticityControl is defined from -45 to 45.  Positive
% ellipticitycontrol indicated CW circular rotation.
% Interactive control is performed with the mouse. Motion with the left
% button down modifies the Tx polarization and right button down modifies
% the Rx polarization.  Horizontal motion changes the angle (Right is
% positive), vertical motion changes the ellipticitycontrol (Up is positive). The
% center button or double-click returns polarization to the default when
% the file was loaded.
%
% INPUTS:
%   filename   - optional : complex image filename(s) that contain a
%                           polarimetric dataset
%                missing  : user is prompted for filename         
%
% OUTPUTS:
%   GUI display, Movie and images as specified
%
% VERSION:
%   1.0
%     - Tim Cox 20100519
%     - initial version, coefficient math provided by Tom Ainsworth (NRL).  MITM tool
%       and readers written by Wade Schwartzkopf (NGA-IDT)
%   1.1
%     - Tim Cox 20100526
%     - Added Analysis Tab w/ Pauli Decomposition and pixel and AOI
%       Angle/Ellip plots
%   1.2
%     - Tim Cox 20100616
%     - Added list box that populates with decomps in "Decomposition"
%       folder
%     - Added three example decomps (Pauli, Simple, DualPolTest)
%     - Added interactive decomposition
%   1.3
%     - Tim Cox 20101202
%     - Added .mat file read (MAT file contains chip data and meta)
%   1.4
%     - Wade Schwartzkopf 20120622
%     - Added embedded mitm_viewer and cleaned up code
%   1.5 
%     - Tim Cox 20200319
%     - Updated to handle RCM (circular transmit polarization) and updated
%       Tab control 
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PolTool_OpeningFcn, ...
                   'gui_OutputFcn',  @PolTool_OutputFcn, ...
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


% --- Executes just before PolTool is made visible.
function PolTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PolTool (see VARARGIN)

% Choose default command line output for PolTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.mitm_hand = hg_mitm_viewer(handles.image);

set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');
set(handles.metaicon,'Color',[0/255 0/255 0/255]);

%format TxPlot
set(handles.TxPlot,'XTick',[]);
set(handles.TxPlot,'YTick',[]);
set(handles.TxPlot,'Box','on');
setAllowAxesZoom(zoom(handles.figure1),handles.TxPlot,false);
setAllowAxesPan(pan(handles.figure1),handles.TxPlot,false);

%format RxPlot
set(handles.RxPlot,'XTick',[]);
set(handles.RxPlot,'YTick',[]);
set(handles.RxPlot,'Box','on');
setAllowAxesZoom(zoom(handles.figure1),handles.RxPlot,false);
setAllowAxesPan(pan(handles.figure1),handles.RxPlot,false);

%set up tab control
%handles.tgroup = uitabgroup('Parent', handles.figure1,'Position', [.017 .015 .4 .30]); When other tabs are visible     
handles.tgroup = uitabgroup('Parent', handles.figure1,'Position', [.017 .015 .72 .30]); %when other tabs are hidden
handles.tab1 = uitab('Parent', handles.tgroup, 'Title', 'Combination');
handles.tab2 = uitab('Parent', handles.tgroup, 'Title', 'Manual Control');
handles.tab3 = uitab('Parent', handles.tgroup, 'Title', 'Save Image/Movie');
handles.tab4 = uitab('Parent', handles.tgroup, 'Title', 'Analysis');

%Place panels into each tab
set(handles.P1,'Parent',handles.tab1)
set(handles.P2,'Parent',handles.tab2)
set(handles.P3,'Parent',handles.tab3)
set(handles.P4,'Parent',handles.tab4)

%Reposition each panel to same location as panel 1
set(handles.P2,'position',get(handles.P1,'position'));
set(handles.P3,'position',get(handles.P1,'position'));
set(handles.P4,'position',get(handles.P1,'position'));

%set default settings
set(handles.AngleControl,'Value',1);
set(handles.EllipticityControl,'Value',1);
set(handles.FrameRate,'String',5);
set(handles.StopMovie,'enable','off');
handles.movieflag = 0;
set(handles.CoAngleEllip,'enable','off');
set(handles.CrossAngleEllip,'enable','off');
set(handles.RxAngleEllip,'enable','off');
set(handles.AngleRes,'String',1);
set(handles.SaveCov,'enable','off');
set(handles.PlotResults,'enable','off');

%transmit polarization.  There are 5 possibilities: HV,H,V,RHC,LHC
handles.TxPol = '';

handles.TxAng = 0;
handles.RxAng = 0;
handles.TxEllip = 0;
handles.RxEllip = 0;

set(handles.TxAngleInc,'String',5);
set(handles.TxEllipInc,'String',5);
set(handles.RxAngleInc,'String',5);
set(handles.RxEllipInc,'String',5);

handles.InMotionFcn = 0;
handles.Select = 'None';
handles.PixelSelect = 0;
handles.SelectingRegion = 0;
handles.PointInImage = 0;
handles.XOffset = 0;
handles.YOffset = 0;
handles.aoi = [];

set(handles.DisplayButtonGroup,'SelectionChangeFcn',@Colormap_Callback);

%populate Decomp Combo box
pathstr=fileparts(mfilename('fullpath'));
filestring1=dir(fullfile(pathstr, 'Decompositions', '*.m'));
filestring2=dir(fullfile(pathstr, 'Decompositions', '*.decomp'));

%trim .m off of decomp functions
for i=1:length(filestring1)
    foo = filestring1(i).name;
    filestring1(i).name = foo(1:length(filestring1(i).name)-2);
end
filestring = cat(1,filestring1,filestring2);
for i=1:length(filestring)
    stringtemp{i} = filestring(i).name;
end
if isempty(filestring)
    stringtemp = 'Pauli';
end
set(handles.DecompCombo,'String',stringtemp);

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[]);
p.addParamValue('segment',1);
p.parse(varargin{:});

if ~isempty(p.Results.filename)
    LoadImage(p.Results.filename,hObject,handles,p.Results.aoi,p.Results.segment);
else
    %set combinations to invisible until files are loaded
    set(handles.Co1,'visible','off');  
    set(handles.Co2,'visible','off');  
    set(handles.Cross1,'visible','off');  
    set(handles.Cross2,'visible','off'); 
    set(handles.CircCo,'visible','off');  
    set(handles.CircCross,'visible','off');  
    set(handles.PauliCo,'visible','off');  
    set(handles.PauliCross,'visible','off'); 
    
    % Update handles structure
    guidata(hObject, handles);
end

% UIWAIT makes PolTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function Colormap_Callback(hObject, eventdata)

handles = guidata(hObject);
if get(handles.Color,'Value')
    colormap(handles.mitm_hand.AxesHandle,jet);
else
    colormap(handles.mitm_hand.AxesHandle,gray);
end

function LoadImage(filenames,hObject,handles,aoi,segment)

% Input arguments
if ischar(filenames) % Input can be single string or cell array of strings
    filenames = {filenames}; % Just treat both types as cell array
end

% Open file(s)
handles.mitm_hand.DataTransformFcn = [];
handles.mitm_hand.close();
if isempty(aoi) % Default is to show entire image
    handles.mitm_hand.openFile(filenames);
else % AOI was passed in
    handles.mitm_hand.openFile(filenames, true);
    pixels_available = floor(getpixelposition(handles.mitm_hand.AxesHandle))-1;
    handles.mitm_hand.setView('CenterPos',aoi(1:2)+(aoi(3:4)/2),...
        'Zoom',max(ceil(aoi(3:4)./pixels_available(3:4))),'Frame',segment(1));
end
meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
MetaIcon(meta,'handle',handles.metaicon);

%determine transmit polarization
pols = cellfun(@(x) split(x,':'), meta.ImageFormation.TxRcvPolarizationProc, 'UniformOutput', false);
TxPol = unique(cellfun(@(x) x{1},pols,'UniformOutput',false));
handles.TxPol = [TxPol{:}];

%draw and label controls based on TxPol (all will have the same Rx Pol)
if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH') %#ok<*BDSCA> %linear Quad
    set(handles.Co1,'String','HH');
    set(handles.Co2,'String','VV');
    set(handles.Cross1,'String','HV');
    set(handles.Cross2,'String','VH');
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0); 
    set(handles.TxAngleMan,'String',0);
    set(handles.RxAngleMan,'String',0);    
    set(handles.TxEllipMan,'String',0);
    set(handles.RxEllipMan,'String',0); 
    set(handles.Co1,'Visible','on');
    set(handles.Cross1,'Visible','on');
    set(handles.Co2,'Visible','on');
    set(handles.Cross2,'Visible','on');
    set(handles.CircCo,'Visible','on');
    set(handles.CircCross,'Visible','on');
    set(handles.PauliCo,'Visible','on');
    set(handles.PauliCross,'Visible','on');
    set(handles.LockPanel,'Visible','on');
    set(handles.TxPolPanel,'Visible','on');
    set(handles.CoAngleEllip,'enable','on');
    set(handles.CrossAngleEllip,'enable','on');
    set(handles.RxAngleEllip,'enable','on');
elseif strcmp(handles.TxPol,'H') %HD
    set(handles.Co1,'String','HH');    
    set(handles.Cross1,'String','HV'); 
    set(handles.Co1,'Visible','on');
    set(handles.Cross1,'Visible','on');
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0); 
    set(handles.Co2,'Visible','off');
    set(handles.Cross2,'Visible','off');
elseif strcmp(handles.TxPol,'V') %VD
    set(handles.Co1,'String','VV');    
    set(handles.Cross1,'String','VH');   
    set(handles.Co1,'Visible','on');
    set(handles.Cross1,'Visible','on');
    set(handles.TxAngle,'String',90);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);  
    handles.TxAng = 90;
    set(handles.Co2,'Visible','off');
    set(handles.Cross2,'Visible','off');
elseif strcmp(handles.TxPol,'RHC') %RHC
    set(handles.Co1,'String','RHC:RHC');    
    set(handles.Cross1,'String','RHC:LHC'); 
    set(handles.Co2,'String','RHC:H');    
    set(handles.Cross2,'String','RHC:V');
    set(handles.Co1,'Visible','on');
    set(handles.Cross1,'Visible','on');
    set(handles.Co2,'Visible','on');
    set(handles.Cross2,'Visible','on');
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',45);
    set(handles.RxEllipticity,'String',45);     
elseif strcmp(handles.TxPol,'LHC') %LHC
    set(handles.Co1,'String','LHC:LHC');    
    set(handles.Cross1,'String','LHC:RHC');
    set(handles.Co2,'String','LHC:H');    
    set(handles.Cross2,'String','LHC:V');
    set(handles.Co1,'Visible','on');
    set(handles.Cross1,'Visible','on');
    set(handles.Co2,'Visible','on');
    set(handles.Cross2,'Visible','on');
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',-45);
    set(handles.RxEllipticity,'String',-45);     
end
if ~strcmp(handles.TxPol,'HV') && ~strcmp(handles.TxPol,'VH')
    set(handles.CircCo,'Visible','off');
    set(handles.CircCross,'Visible','off');
    set(handles.PauliCo,'Visible','off');
    set(handles.PauliCross,'Visible','off');
    set(handles.LockPanel,'Visible','off');
    set(handles.TxPolPanel,'Visible','off');
    set(handles.CoAngleEllip,'enable','off');
    set(handles.CrossAngleEllip,'enable','off');
    set(handles.RxAngleEllip,'enable','on');
    set(handles.RxAngleEllip,'value',1);
end

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

function handles = UpdateImage(handles)

TxAngle = str2double(get(handles.TxAngle,'String'));
TxEllipticity = str2double(get(handles.TxEllipticity,'String'));
RxAngle = str2double(get(handles.RxAngle,'String'));
RxEllipticity = str2double(get(handles.RxEllipticity,'String'));

if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH')
    [coefs(1,1,1),coefs(1,1,2),coefs(1,1,3),coefs(1,1,4)] = ...
    ComputePolCoeff(TxAngle,RxAngle,TxEllipticity,RxEllipticity);
elseif strcmp(handles.TxPol,'H')
    [coefs(1,1,1),coefs(1,1,2)] = ...
    ComputePolCoeff(TxAngle,RxAngle,TxEllipticity,RxEllipticity);
elseif strcmp(handles.TxPol,'V')
    [~,~,coefs(1,1,1),coefs(1,1,2)] = ...
    ComputePolCoeff(TxAngle,RxAngle,TxEllipticity,RxEllipticity);
end

% Apply coefficients
handles.mitm_hand.DataTransformFcn = @(x) sum(bsxfun(@times,x,coefs),3);

if (handles.movieflag == 1)
    %record frame settings
    temp.TxAngle = str2double(get(handles.TxAngle,'String'));
    temp.RxAngle = str2double(get(handles.RxAngle,'String'));
    temp.TxEllip = str2double(get(handles.TxEllipticity,'String'));
    temp.RxEllip = str2double(get(handles.RxEllipticity,'String'));
    temp.Color = get(handles.Color,'Value');
    handles.movieparams = horzcat(handles.movieparams,temp);   
end

% --- Outputs from this function are returned to the command line.
function varargout = PolTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function TxAngle_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TxAngle as text
%        str2double(get(hObject,'String')) returns contents of TxAngle as a double


% --- Executes during object creation, after setting all properties.
function TxAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TxAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TxEllipticity_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipticity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TxEllipticity as text
%        str2double(get(hObject,'String')) returns contents of TxEllipticity as a double


% --- Executes during object creation, after setting all properties.
function TxEllipticity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TxEllipticity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RxAngle_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxAngle as text
%        str2double(get(hObject,'String')) returns contents of RxAngle as a double


% --- Executes during object creation, after setting all properties.
function RxAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RxAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RxEllipticity_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipticity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxEllipticity as text
%        str2double(get(hObject,'String')) returns contents of RxEllipticity as a double


% --- Executes during object creation, after setting all properties.
function RxEllipticity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RxEllipticity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AngleControl.
function AngleControl_Callback(hObject, eventdata, handles)
% hObject    handle to AngleControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AngleControl


% --- Executes on button press in EllipticityControl.
function EllipticityControl_Callback(hObject, eventdata, handles)
% hObject    handle to EllipticityControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EllipticityControl


% --- Executes on button press in Co1.
function Co1_Callback(hObject, eventdata, handles)
% hObject    handle to Co1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH')%#ok<*BDSCA> %linear Quad
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);  
    set(handles.TxAngleMan,'String',0);
    set(handles.RxAngleMan,'String',0);    
    set(handles.TxEllipMan,'String',0);
    set(handles.RxEllipMan,'String',0);  
elseif strcmp(handles.TxPol,'H') %HD   
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);  
    set(handles.TxAngleMan,'String',0);
    set(handles.RxAngleMan,'String',0);    
    set(handles.TxEllipMan,'String',0);
    set(handles.RxEllipMan,'String',0); 
elseif strcmp(handles.TxPol,'V') %VD
     set(handles.TxAngle,'String',90);
    set(handles.RxAngle,'String',90);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);
    set(handles.TxAngleMan,'String',90);
    set(handles.RxAngleMan,'String',90);    
    set(handles.TxEllipMan,'String',0);
    set(handles.RxEllipMan,'String',0); 
    handles.TxAng = 90;
elseif strcmp(handles.TxPol,'RHC') %RHC
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',45);
    set(handles.RxEllipticity,'String',45);  
    set(handles.TxAngleMan,'String',0);
    set(handles.RxAngleMan,'String',0);    
    set(handles.TxEllipMan,'String',45);
    set(handles.RxEllipMan,'String',45); 
elseif strcmp(handles.TxPol,'LHC') %LHC
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',-45);
    set(handles.RxEllipticity,'String',-45);    
    set(handles.TxAngleMan,'String',0);
    set(handles.RxAngleMan,'String',0);    
    set(handles.TxEllipMan,'String',-45);
    set(handles.RxEllipMan,'String',-45); 
end

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Cross1.
function Cross1_Callback(hObject, eventdata, handles)
% hObject    handle to Cross1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH')%#ok<*BDSCA> %linear Quad
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',90);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);   
elseif strcmp(handles.TxPol,'H') %HD   
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',90);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);   
elseif strcmp(handles.TxPol,'V') %VD
     set(handles.TxAngle,'String',90);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);  
    handles.TxAng = 90;
elseif strcmp(handles.TxPol,'RHC') %RHC
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',90);    
    set(handles.TxEllipticity,'String',45);
    set(handles.RxEllipticity,'String',-45);     
elseif strcmp(handles.TxPol,'LHC') %LHC
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',90);    
    set(handles.TxEllipticity,'String',-45);
    set(handles.RxEllipticity,'String',45);     
end

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Co2.
function Co2_Callback(hObject, eventdata, handles)
% hObject    handle to Co2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Quad (VV)
if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH')
    set(handles.TxAngle,'String',90);
    set(handles.RxAngle,'String',90);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0); 
elseif strcmp(handles.TxPol,'RHC')
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',45);
    set(handles.RxEllipticity,'String',0); 
elseif strcmp(handles.TxPol,'LHC')
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',-45);
    set(handles.RxEllipticity,'String',0); 
end

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Cross2.
function Cross2_Callback(hObject, eventdata, handles)
% hObject    handle to Cross2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Quad (VH)
if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH')
    set(handles.TxAngle,'String',90);
    set(handles.RxAngle,'String',00);    
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0); 
elseif strcmp(handles.TxPol,'RHC')
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',45);
    set(handles.RxEllipticity,'String',90); 
elseif strcmp(handles.TxPol,'LHC')
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);    
    set(handles.TxEllipticity,'String',-45);
    set(handles.RxEllipticity,'String',90); 
end

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in CircCo.
function CircCo_Callback(hObject, eventdata, handles)
% hObject    handle to CircCo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%only for Quad
set(handles.TxAngle,'String',0);
set(handles.RxAngle,'String',0);    
set(handles.TxEllipticity,'String',45);
set(handles.RxEllipticity,'String',45); 

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in CircCross.
function CircCross_Callback(hObject, eventdata, handles)
% hObject    handle to CircCross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%only for Quad 
set(handles.TxAngle,'String',0);
set(handles.RxAngle,'String',0);    
set(handles.TxEllipticity,'String',45);
set(handles.RxEllipticity,'String',-45); 

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PauliCo.
function PauliCo_Callback(hObject, eventdata, handles)
% hObject    handle to PauliCo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%only for Quad 
set(handles.TxAngle,'String',135);
set(handles.RxAngle,'String',135);    
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0); 

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PauliCross.
function PauliCross_Callback(hObject, eventdata, handles)
% hObject    handle to PauliCross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%only for Quad 
set(handles.TxAngle,'String',45);
set(handles.RxAngle,'String',135);    
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0); 

PlotPolEllipse(handles.TxPlot,1,handles);
PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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

%get iq filename
[fname, path] = uigetfile( sar_file_extensions('complex'),...
    'Open Image File',pathstr,'MultiSelect', 'on');
if isnumeric(fname)
    return;
end

filenames = {};
if (iscell(fname))
    for i=1:length(fname)
        filenames{i} = strcat(path,fname{i});
    end
else      
    filenames{1} = strcat(path,fname);
end

setpref('matlab_sar_toolbox','last_used_directory',path); %store path

%Load Image
LoadImage(filenames,hObject,handles,[],1);



function PlotPolEllipse(handle,Tx,handles)

%Tx flag 1 for Tx, 0 for Rx
if Tx
    Ang = str2double(get(handles.TxAngle,'String'));
    Ellip = str2double(get(handles.TxEllipticity,'String'));
else
    Ang = str2double(get(handles.RxAngle,'String'));
    Ellip = str2double(get(handles.RxEllipticity,'String'));
end

Ang = -1*Ang;

%first thing to do is delete any plots that exist within the givin handle
delete(get(handle,'Children'));

%we'll make the radius of this circle 10
%compute the semi-major/semi-minor angle
SemiMajor = 10*cosd(Ellip);
SemiMinor = 10*sind(abs(Ellip));

% Points for Ellipse
t=linspace(0,2*pi,100);
x=SemiMajor*cos(t);
y=SemiMinor*sin(t);
[th,r]=cart2pol(x,y);
[x,y]=pol2cart(th-Ang*pi/180,r);

% Plot Ellipse
axes(handle);
hold on;
plot(x,y,'LineWidth',3,'Color',[0 0.4470 0.7410]); % Blue
xlim([-10 10]);
ylim([-10 10]);
set(handle,'XTick',[]);
set(handle,'YTick',[]);

%plot vertical and horizontal lines for reference
plot([-10 10],[0 0],'Color','k','LineWidth',1);
plot([0 0],[-10 10],'Color','k','LineWidth',1);

%Draw Rotation Angle based on sign of Ellipticity. If Ellipticity is small
%then we'll draw a two headed arrow rotated for the specified angle. We'll
%draw this in the lower left or right corner, depending on the angle.

if (Ang < 0)
    Ang = Ang + 360;
end

if (abs(Ellip) > 2.5)
    %draw rotation direction indicator
    if (Ellip > 0)
        %circle part
        count = 0;
        for i=0:5:270
            count = count + 1;
            circX(count) = 2*cosd(i);
            circY(count) = 2*sind(i);
        end
        %Arrows
        Arrow1X = [2 2.75]; 
        Arrow1Y = [0 0.75];
        Arrow2X = [2 1.25]; 
        Arrow2Y = [0 0.75];        
    else
        %circle part
        count = 0;
        for i=-90:5:180
            count = count + 1;
            circX(count) = 2*cosd(i);
            circY(count) = 2*sind(i);
        end
        %Arrows
        Arrow1X = [-2 -2.75]; 
        Arrow1Y = [0 0.75];
        Arrow2X = [-2 -1.25]; 
        Arrow2Y = [0 0.75];  
    end
    %move based on Angle
    if ((Ang > 0 && Ang < 90) || ...
        (Ang > 180 && Ang < 270))
        %draw rotation indicator in lower left corner
        circX = circX - 7.5; 
        circY = circY - 7.5; 
        Arrow1X = Arrow1X -7.5;
        Arrow1Y = Arrow1Y -7.5;
        Arrow2X = Arrow2X -7.5;
        Arrow2Y = Arrow2Y -7.5;
    else
        %draw rotation indicator in lower right corner
        circX = circX + 7.5; 
        circY = circY - 7.5; 
        Arrow1X = Arrow1X +7.5;
        Arrow1Y = Arrow1Y -7.5;
        Arrow2X = Arrow2X +7.5;
        Arrow2Y = Arrow2Y -7.5;
    end

    plot(circX,circY,'Color','k','LineWidth',2);
    plot(Arrow1X,Arrow1Y,'Color','k','LineWidth',2);
    plot(Arrow2X,Arrow2Y,'Color','k','LineWidth',2);
else
    %draw two headed arrow rotated to angle
    %Arrow will consist of 5 lines
    Line1X = [0 0];
    Line1Y = [-2 2];
    Line2X = [0 -0.75];
    Line2Y = [2 1.25];
    Line3X = [0 0.75];
    Line3Y = [2 1.25]; 
    Line4X = [0 -0.75];
    Line4Y = [-2 -1.25];
    Line5X = [0 0.75];
    Line5Y = [-2 -1.25];
        
    %rotate lines
    [Line1Th Line1R] = cart2pol(Line1X,Line1Y);
    [Line1X Line1Y] = pol2cart(Line1Th - Ang*pi/180 + pi/2,Line1R);
    [Line2Th Line2R] = cart2pol(Line2X,Line2Y);
    [Line2X Line2Y] = pol2cart(Line2Th - Ang*pi/180 + pi/2,Line2R);
    [Line3Th Line3R] = cart2pol(Line3X,Line3Y);
    [Line3X Line3Y] = pol2cart(Line3Th - Ang*pi/180 + pi/2,Line3R);
    [Line4Th Line4R] = cart2pol(Line4X,Line4Y);
    [Line4X Line4Y] = pol2cart(Line4Th - Ang*pi/180 + pi/2,Line4R);
    [Line5Th Line5R] = cart2pol(Line5X,Line5Y);
    [Line5X Line5Y] = pol2cart(Line5Th - Ang*pi/180 + pi/2,Line5R);
    
    
    %move based on Angle
    if ((Ang > 0 && Ang < 90) || ...
        (Ang > 180 && Ang < 270))
        Line1X = Line1X -7.5;
        Line1Y = Line1Y -7.5;
        Line2X = Line2X -7.5;
        Line2Y = Line2Y -7.5;
        Line3X = Line3X -7.5;
        Line3Y = Line3Y -7.5;
        Line4X = Line4X -7.5;
        Line4Y = Line4Y -7.5;
        Line5X = Line5X -7.5;
        Line5Y = Line5Y -7.5;
    else
        Line1X = Line1X +7.5;
        Line1Y = Line1Y -7.5;
        Line2X = Line2X +7.5;
        Line2Y = Line2Y -7.5;
        Line3X = Line3X +7.5;
        Line3Y = Line3Y -7.5;
        Line4X = Line4X +7.5;
        Line4Y = Line4Y -7.5;
        Line5X = Line5X +7.5;
        Line5Y = Line5Y -7.5;
    end
        
    plot(Line1X,Line1Y,'Color','k','LineWidth',2);
    plot(Line2X,Line2Y,'Color','k','LineWidth',2);
    plot(Line3X,Line3Y,'Color','k','LineWidth',2);
    plot(Line4X,Line4Y,'Color','k','LineWidth',2);
    plot(Line5X,Line5Y,'Color','k','LineWidth',2);
end

hold off;

%make sure image is updated now
drawnow expose;


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isfield(handles,'Select')&&strcmp(handles.Select,'None'))
    return;
end

if ~isfield(handles,'InMotionFcn') || handles.InMotionFcn == 1 || ...
        isempty(handles.mitm_hand.Metadata) || ...
        handles.SelectingRegion == 1 || handles.PixelSelect == 1
    return;
end

handles.InMotionFcn = 1;
% Update handles structure
guidata(hObject, handles);

AngleControl = get(handles.AngleControl,'Value');
EllipControl = get(handles.EllipticityControl,'Value');

AngleLock = get(handles.AngleLock,'Value');
EllipLock = get(handles.EllipticityLock,'Value');

%Transmit
TxAng = str2double(get(handles.TxAngle,'String')); 
TxEllip = str2double(get(handles.TxEllipticity,'String')); 
%Recieve
RxAng = str2double(get(handles.RxAngle,'String')); 
RxEllip = str2double(get(handles.RxEllipticity,'String')); 

pos = get(gcf,'CurrentPoint');
delta = handles.pos - pos;
if (AngleControl)
    DeltaAngle = -1*delta(1)*2;
else
    DeltaAngle = 0;
end
if (EllipControl)
    DeltaEllip = delta(2)*2;   
else
    DeltaEllip = 0;
end

if (sqrt(DeltaAngle*DeltaAngle+DeltaEllip*DeltaEllip) < 1)
    handles.InMotionFcn = 0;
    % Update handles structure
    guidata(hObject, handles);
    return;
end

if (strcmp(handles.Select,'Tx'))
    TxAng = TxAng - DeltaAngle;
    if (AngleLock)
        RxAng = RxAng - DeltaAngle;
    end
    TxEllip = TxEllip - DeltaEllip;
    if (EllipLock)
        RxEllip = RxEllip - DeltaEllip;
    end
else
    RxAng = RxAng - DeltaAngle;
    if (AngleLock)
        TxAng = TxAng - DeltaAngle;
    end
    RxEllip = RxEllip - DeltaEllip;
    if (EllipLock)
        TxEllip = TxEllip - DeltaEllip;
    end
end

%Angle is between 0 and 360 
%Ellipticity is between -45 and 45
while (TxAng > 360)
    TxAng = TxAng - 360;
end    
while (TxAng < 0)
    TxAng = TxAng + 360;
end    
if (TxEllip > 45)
    TxEllip = TxEllip - 90;
end    
if (TxEllip < -45)
    TxEllip = TxEllip + 90;
end
while (RxAng > 360)
    RxAng = RxAng - 360;
end    
while (RxAng < 0)
    RxAng = RxAng + 360;
end    
if (RxEllip > 45)
    RxEllip = RxEllip - 90;
end    
if (RxEllip < -45)
    RxEllip = RxEllip + 90;
end
    
set(handles.TxAngle,'String',round(TxAng));
set(handles.TxEllipticity,'String',round(TxEllip));
set(handles.RxAngle,'String',round(RxAng));
set(handles.RxEllipticity,'String',round(RxEllip));

handles.pos = pos;

%update Ellipse Plots (only if necessary)
if (strcmp(handles.Select,'Tx') || AngleLock || EllipLock)
    PlotPolEllipse(handles.TxPlot,1,handles);
end
if (strcmp(handles.Select,'Rx') || AngleLock || EllipLock)
    PlotPolEllipse(handles.RxPlot,0,handles);
end

%Update Image
handles = UpdateImage(handles);

handles.InMotionFcn = 0;

% Update handles structure
guidata(hObject, handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

typ = get( gcf, 'SelectionType' );

if strcmpi(typ,'normal')
    if strcmpi(handles.TxPol,'HV')
        %Left Mouse Click
        handles.Select = 'Tx';
    else
        handles.Select = 'None'; %dual pol
    end
elseif strcmpi(typ,'alt')
    %right mouse click
    handles.Select = 'Rx';
else
    %double click, reset to default view
    Co1_Callback(hObject, eventdata, handles)
end

%store initial point
handles.pos = get(gcf,'CurrentPoint'); 

% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Select = 'None';

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in EllipticityLock.
function EllipticityLock_Callback(hObject, eventdata, handles)
% hObject    handle to EllipticityLock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EllipticityLock


% --- Executes on button press in AngleLock.
function AngleLock_Callback(hObject, eventdata, handles)
% hObject    handle to AngleLock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AngleLock



function TxAngleMan_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngleMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TxAngleMan as text
%        str2double(get(hObject,'String')) returns contents of TxAngleMan as a double


% --- Executes during object creation, after setting all properties.
function TxAngleMan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TxAngleMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TxAngleInc_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngleInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TxAngleInc as text
%        str2double(get(hObject,'String')) returns contents of TxAngleInc as a double


% --- Executes during object creation, after setting all properties.
function TxAngleInc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TxAngleInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TxAngleUp.
function TxAngleUp_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngleUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxAngle = str2double(get(handles.TxAngleMan,'String'));
TxInc = str2double(get(handles.TxAngleInc,'String'));

TxAngle = TxAngle+TxInc;
if TxAngle > 360; TxAngle = TxAngle-360; end

set(handles.TxAngleMan,'String',round(TxAngle));
set(handles.TxAngle,'String',round(TxAngle));

PlotPolEllipse(handles.TxPlot,1,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in TxAngleDown.
function TxAngleDown_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngleDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxAngle = str2double(get(handles.TxAngleMan,'String'));
TxInc = str2double(get(handles.TxAngleInc,'String'));

TxAngle = TxAngle-TxInc;
if TxAngle < 0; TxAngle = TxAngle+360; end

set(handles.TxAngleMan,'String',round(TxAngle));
set(handles.TxAngle,'String',round(TxAngle));

PlotPolEllipse(handles.TxPlot,1,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

function TxEllipMan_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TxEllipMan as text
%        str2double(get(hObject,'String')) returns contents of TxEllipMan as a double


% --- Executes during object creation, after setting all properties.
function TxEllipMan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TxEllipMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TxEllipInc_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TxEllipInc as text
%        str2double(get(hObject,'String')) returns contents of TxEllipInc as a double


% --- Executes during object creation, after setting all properties.
function TxEllipInc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TxEllipInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TxEllipUp.
function TxEllipUp_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxEllip = str2double(get(handles.TxEllipMan,'String'));
TxInc = str2double(get(handles.TxEllipInc,'String'));

TxEllip = TxEllip+TxInc;
if TxEllip > 45; TxEllip = TxEllip-90; end

set(handles.TxEllipMan,'String',round(TxEllip));
set(handles.TxEllipticity,'String',round(TxEllip));

PlotPolEllipse(handles.TxPlot,1,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in TxEllipDown.
function TxEllipDown_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxEllip = str2double(get(handles.TxEllipMan,'String'));
TxInc = str2double(get(handles.TxEllipInc,'String'));

TxEllip = TxEllip-TxInc;
if TxEllip < -45; TxEllip = TxEllip+90; end

set(handles.TxEllipMan,'String',round(TxEllip));
set(handles.TxEllipticity,'String',round(TxEllip));

PlotPolEllipse(handles.TxPlot,1,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in RxEllipDown.
function RxEllipDown_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxEllip = str2double(get(handles.RxEllipMan,'String'));
RxInc = str2double(get(handles.RxEllipInc,'String'));

RxEllip = RxEllip-RxInc;
if RxEllip < -45; RxEllip = RxEllip+90; end

set(handles.RxEllipMan,'String',round(RxEllip));
set(handles.RxEllipticity,'String',round(RxEllip));

PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in RxEllipUp.
function RxEllipUp_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxEllip = str2double(get(handles.RxEllipMan,'String'));
RxInc = str2double(get(handles.RxEllipInc,'String'));

RxEllip = RxEllip+RxInc;
if RxEllip > 45; RxEllip = RxEllip-90; end

set(handles.RxEllipMan,'String',round(RxEllip));
set(handles.RxEllipticity,'String',round(RxEllip));

PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);


function RxEllipInc_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxEllipInc as text
%        str2double(get(hObject,'String')) returns contents of RxEllipInc as a double


% --- Executes during object creation, after setting all properties.
function RxEllipInc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RxEllipInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RxEllipMan_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxEllipMan as text
%        str2double(get(hObject,'String')) returns contents of RxEllipMan as a double


% --- Executes during object creation, after setting all properties.
function RxEllipMan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RxEllipMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RxAngleDown.
function RxAngleDown_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngleDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxAngle = str2double(get(handles.RxAngleMan,'String'));
RxInc = str2double(get(handles.RxAngleInc,'String'));

RxAngle = RxAngle-RxInc;
if RxAngle < 0; RxAngle = RxAngle+360; end

set(handles.RxAngleMan,'String',round(RxAngle));
set(handles.RxAngle,'String',round(RxAngle));

PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in RxAngleUp.
function RxAngleUp_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngleUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxAngle = str2double(get(handles.RxAngleMan,'String'));
RxInc = str2double(get(handles.RxAngleInc,'String'));

RxAngle = RxAngle+RxInc;
if RxAngle > 360; RxAngle = RxAngle-360; end

set(handles.RxAngleMan,'String',round(RxAngle));
set(handles.RxAngle,'String',round(RxAngle));

PlotPolEllipse(handles.RxPlot,0,handles);

handles = UpdateImage(handles);

% Update handles structure
guidata(hObject, handles);


function RxAngleInc_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngleInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxAngleInc as text
%        str2double(get(hObject,'String')) returns contents of RxAngleInc as a double


% --- Executes during object creation, after setting all properties.
function RxAngleInc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RxAngleInc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RxAngleMan_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngleMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxAngleMan as text
%        str2double(get(hObject,'String')) returns contents of RxAngleMan as a double


% --- Executes during object creation, after setting all properties.
function RxAngleMan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RxAngleMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MovieName_Callback(hObject, eventdata, handles)
% hObject    handle to MovieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MovieName as text
%        str2double(get(hObject,'String')) returns contents of MovieName as a double


% --- Executes during object creation, after setting all properties.
function MovieName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MovieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseMovie.
function BrowseMovie_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseMovie (see GCBO)
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

[fname, path] = uiputfile( {'*.gif','GIF Files (*.gif)'},'Save Animated GIF',pathstr);

setpref('matlab_sar_toolbox','last_used_directory',path); %store path

set(handles.MovieName,'String',strcat(path,fname));
set(handles.StartMovie,'Enable','on');


% --- Executes on button press in SnapImage.
function SnapImage_Callback(hObject, eventdata, handles)
% hObject    handle to SnapImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

F = getframe(handles.mitm_hand.AxesHandle );
dirname = fileparts(handles.mitm_hand.Filename{1});
outfile = uiputfile('*.jpg','Save Image Snap',dirname);
imwrite(F.cdata,[dirname filesep outfile],'jpg');   
msgbox('Image Saved');

% --- Executes on button press in StartMovie.
function StartMovie_Callback(hObject, eventdata, handles)
% hObject    handle to StartMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%we'll just store the Pol values and zoom/contrast/color params and then
%create the movie 
handles.movieparams = [];

%set avi recording flag
handles.movieflag = 1;

set(handles.StartMovie,'Enable','off');
set(handles.StopMovie,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in StopMovie.
function StopMovie_Callback(hObject, eventdata, handles)
% hObject    handle to StopMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.movieflag = 0;

set(handles.StartMovie,'Enable','off');
set(handles.StopMovie,'Enable','off');

framerate = str2double(get(handles.FrameRate,'String'));
moviename = get(handles.MovieName,'String');

numframes = length(handles.movieparams);
set(handles.TotFrames,'String',numframes);

for i=1:numframes
    set(handles.TxAngle,'String',round(handles.movieparams(i).TxAngle));
    set(handles.RxAngle,'String',round(handles.movieparams(i).RxAngle));
    set(handles.TxEllipticity,'String',round(handles.movieparams(i).TxEllip));
    set(handles.RxEllipticity,'String',round(handles.movieparams(i).RxEllip));
    set(handles.Color,'Value',handles.movieparams(i).Color);
    UpdateImage(handles);
    PlotPolEllipse(handles.RxPlot,0,handles);
    PlotPolEllipse(handles.TxPlot,1,handles);
    F = getframe(handles.mitm_hand.AxesHandle);
    if i==1
        [X,map] = rgb2ind(F.cdata,256);
        imwrite(X,map,moviename,'gif','LoopCount',Inf,'DelayTime',1/framerate);
    else
        X = rgb2ind(F.cdata,map);
        imwrite(X,map,moviename,'gif','WriteMode','append','DelayTime',1/framerate);
    end
       
    %update status on GUI
    set(handles.CurFrame,'String',i);
end

set(handles.TotFrames,'String','');
set(handles.CurFrame,'String','');

set(handles.StartMovie,'Enable','on');
set(handles.StopMovie,'Enable','off');

% Update handles structure
guidata(hObject, handles);

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


% --- Executes on selection change in DecompCombo.
function DecompCombo_Callback(hObject, eventdata, handles)
% hObject    handle to DecompCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DecompCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DecompCombo


% --- Executes during object creation, after setting all properties.
function DecompCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DecompCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunDecomp.
function RunDecomp_Callback(hObject, eventdata, handles)
% hObject    handle to RunDecomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%put in zeros for any pols that have no data associated with them
data_size = size(handles.mitm_hand.ComplexData);
data_size = data_size([1 2]);

if strcmp(handles.TxPol,'HV') || strcmp(handles.TxPol,'VH')
    HHData = single(handles.mitm_hand.ComplexData(:,:,handles.HH));
    HVData = single(handles.mitm_hand.ComplexData(:,:,handles.HV));
    VHData = single(handles.mitm_hand.ComplexData(:,:,handles.VH));
    VVData = single(handles.mitm_hand.ComplexData(:,:,handles.VV));
elseif strcmp(handles.TxPol,'H')
    HHData = single(handles.mitm_hand.ComplexData(:,:,handles.HH));
    HVData = single(handles.mitm_hand.ComplexData(:,:,handles.HV));
    VHData = zeros(data_size);
    VVData = zeros(data_size);
elseif strcmp(handles.TxPol,'V')
    VHData = single(handles.mitm_hand.ComplexData(:,:,handles.VH));
    VVData = single(handles.mitm_hand.ComplexData(:,:,handles.VV));
    HHData = zeros(data_size);
    HVData = zeros(data_size);
elseif strcmp(handles.TxPol,'RHC') || strcmp(handles.TxPol,'LHC')
    VHData = single(handles.mitm_hand.ComplexData(:,:,handles.VH));
    VVData = single(handles.mitm_hand.ComplexData(:,:,handles.VV));
    HHData = zeros(data_size);
    HVData = zeros(data_size);
end

%run decomposition from list
%get name of function to call
decompindex = get(handles.DecompCombo,'Value');
decompstrings = get(handles.DecompCombo,'String');
decomp = decompstrings(decompindex);

%determine if this is a matlab file or a decomp file (decomp file has
%extension)
if (isempty(strfind(decomp{1},'.decomp')))
    RGBData = feval(decomp{1}, HHData, HVData, VHData, VVData);
else
    pathstr = fileparts(which('PolTool'));
    DecompDir = [pathstr filesep 'Decompositions' filesep];

    RGBData = CreateDecomp([DecompDir decomp{1}], HHData, HVData, VHData, VVData);
end

if isdeployed
    %imtool does not work in compiled Matlab
    figure;
    imagesc(RGBData);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
else
    imtool(RGBData);
end



function RedString_Callback(hObject, eventdata, handles)
% hObject    handle to RedString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RedString as text
%        str2double(get(hObject,'String')) returns contents of RedString as a double


% --- Executes during object creation, after setting all properties.
function RedString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetRed.
function SetRed_Callback(hObject, eventdata, handles)
% hObject    handle to SetRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.RedTxAngle = str2double(get(handles.TxAngle,'String'));
handles.RedTxEllip = str2double(get(handles.TxEllipticity,'String'));
handles.RedRxAngle = str2double(get(handles.RxAngle,'String'));
handles.RedRxEllip = str2double(get(handles.RxEllipticity,'String'));

PolString = sprintf('Tx: %d/%d, Rx: %d/%d',handles.RedTxAngle,...
    handles.RedTxEllip,handles.RedRxAngle,handles.RedRxEllip);
set(handles.RedString,'String',PolString);
guidata(hObject, handles);

function GreenString_Callback(hObject, eventdata, handles)
% hObject    handle to GreenString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GreenString as text
%        str2double(get(hObject,'String')) returns contents of GreenString as a double


% --- Executes during object creation, after setting all properties.
function GreenString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GreenString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetGreen.
function SetGreen_Callback(hObject, eventdata, handles)
% hObject    handle to SetGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.GreenTxAngle = str2double(get(handles.TxAngle,'String'));
handles.GreenTxEllip = str2double(get(handles.TxEllipticity,'String'));
handles.GreenRxAngle = str2double(get(handles.RxAngle,'String'));
handles.GreenRxEllip = str2double(get(handles.RxEllipticity,'String'));

PolString = sprintf('Tx: %d/%d, Rx: %d/%d',handles.GreenTxAngle,...
    handles.GreenTxEllip,handles.GreenRxAngle,handles.GreenRxEllip);
set(handles.GreenString,'String',PolString);
guidata(hObject, handles);

function BlueString_Callback(hObject, eventdata, handles)
% hObject    handle to BlueString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BlueString as text
%        str2double(get(hObject,'String')) returns contents of BlueString as a double


% --- Executes during object creation, after setting all properties.
function BlueString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlueString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SetBlue.
function SetBlue_Callback(hObject, eventdata, handles)
% hObject    handle to SetBlue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.BlueTxAngle = str2double(get(handles.TxAngle,'String'));
handles.BlueTxEllip = str2double(get(handles.TxEllipticity,'String'));
handles.BlueRxAngle = str2double(get(handles.RxAngle,'String'));
handles.BlueRxEllip = str2double(get(handles.RxEllipticity,'String'));

PolString = sprintf('Tx: %d/%d, Rx: %d/%d',handles.BlueTxAngle,...
    handles.BlueTxEllip,handles.BlueRxAngle,handles.BlueRxEllip);
set(handles.BlueString,'String',PolString);
guidata(hObject, handles);


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%generate "custom" defined decomposition

data = handles.mitm_hand.ComplexData;

%Red Channel
[coefs(1,1,handles.HH) coefs(1,1,handles.HV) ...
    coefs(1,1,handles.VH) coefs(1,1,handles.VV)] = ...
    ComputePolCoeff(handles.RedTxAngle,handles.RedRxAngle,...
    handles.RedTxEllip,handles.RedRxEllip);
% Apply coefficients
RedData = abs(sum(bsxfun(@times,data,coefs),3));

%Green Channel
[coefs(1,1,handles.HH) coefs(1,1,handles.HV) ...
    coefs(1,1,handles.VH) coefs(1,1,handles.VV)] = ...
    ComputePolCoeff(handles.GreenTxAngle,handles.GreenRxAngle,...
    handles.GreenTxEllip,handles.GreenRxEllip);
% Apply coefficients
GreenData = abs(sum(bsxfun(@times,data,coefs),3));

%Blue Channel
[coefs(1,1,handles.HH) coefs(1,1,handles.HV) ...
    coefs(1,1,handles.VH) coefs(1,1,handles.VV)] = ...
    ComputePolCoeff(handles.BlueTxAngle,handles.BlueRxAngle,...
    handles.BlueTxEllip,handles.BlueRxEllip);
% Apply coefficients
BlueData = abs(sum(bsxfun(@times,data,coefs),3));

%scale data to mean + 3*sigma 
RedCutoff = mean(RedData(:)) + 3*std(RedData(:));
GreenCutoff = mean(GreenData(:)) + 3*std(GreenData(:));
BlueCutoff = mean(BlueData(:)) + 3*std(BlueData(:));

RedData(RedData > RedCutoff) = RedCutoff;
GreenData(GreenData > GreenCutoff) = GreenCutoff;
BlueData(BlueData > BlueCutoff) = BlueCutoff;

RedData = RedData./RedCutoff;
GreenData = GreenData./GreenCutoff;
BlueData = BlueData./BlueCutoff;

RGBData(:,:,1) = RedData;
RGBData(:,:,2) = GreenData;
RGBData(:,:,3) = BlueData;

if isdeployed
    %imtool does not work in compiled Matlab
    figure;
    imagesc(RGBData);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
else
    imtool(RGBData);
end



% --- Executes on button press in SaveDecomp.
function SaveDecomp_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDecomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[pathstr] = fileparts(which('PolTool'));
DecompDir = sprintf('%s/Decompositions',pathstr);

%get name of decomp
[name pathstr] = uiputfile('*.decomp','Save Decomposition Settings',DecompDir);

DecompName = strcat(pathstr,name);

fid = fopen(DecompName,'w');

fprintf(fid,'RedTxAngle = %d\n',handles.RedTxAngle);
fprintf(fid,'RedTxEllip = %d\n',handles.RedTxEllip);
fprintf(fid,'RedRxAngle = %d\n',handles.RedRxAngle);
fprintf(fid,'RedRxEllip = %d\n\n',handles.RedRxEllip);

fprintf(fid,'GreenTxAngle = %d\n',handles.GreenTxAngle);
fprintf(fid,'GreenTxEllip = %d\n',handles.GreenTxEllip);
fprintf(fid,'GreenRxAngle = %d\n',handles.GreenRxAngle);
fprintf(fid,'GreenRxEllip = %d\n\n',handles.GreenRxEllip);

fprintf(fid,'BlueTxAngle = %d\n',handles.BlueTxAngle);
fprintf(fid,'BlueTxEllip = %d\n',handles.BlueTxEllip);
fprintf(fid,'BlueRxAngle = %d\n',handles.BlueRxAngle);
fprintf(fid,'BlueRxEllip = %d\n',handles.BlueRxEllip);

fclose(fid);

%update list box with new decomp file
%populate decomp list (from Decompositions folder)
[pathstr] = fileparts(which('PolTool'));
SearchDir = sprintf('%s/Decompositions/*.m',pathstr);
files = dir(SearchDir);
NumMFiles = length(files);
for i=1:NumMFiles
    temp = files(i).name;
    filestring{i} = temp(1:length(temp)-2);
end

SearchDir = sprintf('%s/Decompositions/*.decomp',pathstr);
files = dir(SearchDir);
NumDecompFiles = length(files);
for i=(NumMFiles+1):(NumMFiles+NumDecompFiles)
    filestring{i} = files(i-NumMFiles).name;    
end

set(handles.DecompCombo,'String',filestring);


% --- Executes on button press in SelectPixel.
function SelectPixel_Callback(hObject, eventdata, handles)
% hObject    handle to SelectPixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%determine what mode we're in
label = get(handles.SelectPixel,'String');

if strcmpi(label,'Select Pixel')
    set(handles.SelectPixel,'String','Delete Point');
    set(handles.SelectAOI,'Enable','off');
    handles.PixelSelect = 1;
    guidata(hObject, handles);
    %place impoint at center of image
    r = getpixelposition(handles.mitm_hand.AxesHandle);
    handles.point = impoint(handles.mitm_hand.AxesHandle,round(r(3)/2),round(r(4)/2));
    kids = get(handles.point,'Children');
    set(kids(1),'MarkerSize',15);
    set(kids(1),'MarkerEdgeColor','r');
    set(kids(1),'LineWidth',2)
    set(kids(2),'MarkerFaceColor','y')
    set(kids(2),'MarkerSize',12)
    set(kids(2),'LineWidth',2)
    handles.PointInImage = 1;
    set(handles.PlotResults,'Enable','on');
    set(handles.SaveCov,'Enable','on');
    guidata(hObject, handles);
else
    handles.PointInImage = 0;
    handles.PixelSelect = 0;
    delete(handles.point);
    set(handles.SelectPixel,'String','Select Pixel');
    set(handles.SelectAOI,'Enable','on');
    set(handles.PlotResults,'Enable','off');
    set(handles.SaveCov,'Enable','off');
    guidata(hObject, handles);
end

% --- Executes on button press in SelectAOI.
function SelectAOI_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

label = get(handles.SelectAOI,'String');

if strcmpi(label,'Select AOI')
    set(handles.SelectAOI,'String','Delete AOI');
    set(handles.SelectPixel,'Enable','off');
    handles.PixelSelect = 1;
    guidata(hObject, handles);
    r = getpixelposition(handles.mitm_hand.AxesHandle);
    pos = [r(3)/4 r(4)/4 r(3)/2-1 r(4)/2-1];
    handles.rect = imrect(handles.mitm_hand.AxesHandle,pos);
    setColor(handles.rect,'y');
    set(handles.PlotResults,'Enable','on');
    set(handles.SaveCov,'Enable','on');
    guidata(hObject, handles);
else
    handles.PixelSelect = 0;
    delete(handles.rect);
    set(handles.SelectAOI,'String','Select AOI');
    set(handles.SelectPixel,'Enable','on');
    set(handles.PlotResults,'Enable','off');
    set(handles.SaveCov,'Enable','off');
    guidata(hObject, handles);
end

% --- Executes on button press in SaveCov.
function SaveCov_Callback(hObject, eventdata, handles)
% hObject    handle to SaveCov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cov_mat = get_aoi_covariance(handles);

%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

[fname, path] = uiputfile( {'*.cov','COV Files (*.cov)'},'Save Covariance File',pathstr);
if isnumeric(fname), return, end; % Cancel

%write data to file
fid = fopen(strcat(path,fname),'w');
pol_order = handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageFormation.TxRcvPolarizationProc;
fprintf(fid,'Covariance: %s',pol_order{1});
for i=2:numel(pol_order)
    fprintf(fid,',%s',pol_order{i});
end
fprintf(fid,'\n');
for i=1:size(cov_mat,1)
    for j=1:size(cov_mat,2)
        fprintf(fid,'Cov(%d,%d): %21.18f %21.18f\n',i,j,real(cov_mat(i,j)),imag(cov_mat(i,j)));
    end
end
fclose(fid);

% --- Executes on button press in PlotResults.
function PlotResults_Callback(hObject, eventdata, handles)
% hObject    handle to PlotResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AngRes = str2double(get(handles.AngleRes,'String'));
ylabels = -45:AngRes:45;
xlabels = 0:AngRes:180;
if get(handles.CoAngleEllip,'Value') % Co-Pol Angle/Ellipse Plot
    TxAngle = xlabels;
    RxAngle = xlabels;
    TxEllip = ylabels;
    RxEllip = ylabels;
    title_str = 'Co-Pol Angle/Ellipticity Plot';
elseif (get(handles.CrossAngleEllip,'Value')) % Cross-Pol Angle/Ellipse Plot
    TxAngle = xlabels;
    RxAngle = xlabels+90;
    TxEllip = ylabels;
    RxEllip = -ylabels;
    title_str = 'Cross-Pol Angle/Ellipticity Plot';
else % Rx - Angle/Ellip
    TxAngle = ones(1,numel(xlabels))*str2double(get(handles.TxAngle,'String'));
    TxEllip = ones(1,numel(ylabels))*str2double(get(handles.TxEllipticity,'String'));
    RxAngle = xlabels;
    RxEllip = ylabels;
    title_str = 'Rx Angle/Ellipticity Plot';
end
coefs = zeros(1,size(handles.mitm_hand.ComplexData,3));
data = zeros(numel(ylabels),numel(xlabels));
cov_mat = get_aoi_covariance(handles);
for i = 1:numel(ylabels)
    for j = 1:numel(xlabels)
        [coefs(handles.HH) coefs(handles.HV) ...
            coefs(handles.VH) coefs(handles.VV)] = ...
            ComputePolCoeff(TxAngle(j),RxAngle(j),TxEllip(i),RxEllip(i));
        data(i,j) = abs(coefs*cov_mat*coefs');
    end
end
figure;
imagesc(xlabels,ylabels,data);colormap(jet);
xlabel('Angle (deg)');
ylabel('Ellipticity (deg)');
xlim(xlabels([1 end]));
ylim(ylabels([1 end]));
title(title_str);

guidata(hObject, handles);

% Get covariance from selected point or AOI
function cov_mat = get_aoi_covariance(handles)
if (handles.PointInImage == 1)
    % Get position of point
    pos = round(getPosition(handles.point));
    % Create psuedo-covariance matrix for a single pixel.  Obviously a
    % single value has no variance, but we assume a zero-mean for SAR data.
    complex_data = single(squeeze(handles.mitm_hand.ComplexData(pos(2),pos(1),:)))';
else % AOI
    screen_coords = handles.rect.getPosition();
    file_coords = handles.mitm_hand.axescoords2native(screen_coords(1:2));
    file_coords(3:4) = screen_coords(3:4)/handles.mitm_hand.Zoom;
    xmin = max(1,round(file_coords(1)));
    xmax = round(file_coords(1)+file_coords(3));
    ymin = max(1,round(file_coords(2)));
    ymax = round(file_coords(2)+file_coords(4));
    complex_data = handles.mitm_hand.readerobj{handles.mitm_hand.Frame}.read_chip...
        ([xmin xmax],[ymin ymax]);
    complex_data = single(complex_data); % Required for doing math on data
    % Reshape data so that each column is a variable (polarimetric channel)
    % and each row is an observation (pixel).
    complex_data = reshape(complex_data,size(complex_data,1)*size(complex_data,2),size(complex_data,3));
end
% We do not demean here as some covariance computations would (like
% MATLAB's COV function) since SAR data is assumed to be zero-mean.
cov_mat = complex_data'*complex_data;


function AngleRes_Callback(hObject, eventdata, handles)
% hObject    handle to AngleRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AngleRes as text
%        str2double(get(hObject,'String')) returns contents of AngleRes as a double


% --- Executes during object creation, after setting all properties.
function AngleRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AngleRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
