function varargout = Taser(varargin)
% Taser Main routine to launch interactive matlab tools
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Main entry point to SAR Matlab Toolbox.  Allows user to select file and 
% AOI and then launch desired analysis algorithm.  Algorithms are populated
% by an external text file (Algorithms.txt).  

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Taser_OpeningFcn, ...
                   'gui_OutputFcn',  @Taser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    warning off all;
    gui_State.gui_Callback = str2func(varargin{1});
    warning on all;
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Taser is made visible.
function Taser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Taser (see VARARGIN)


% Choose default command line output for Taser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.mitm_hand = hg_mitm_viewer(handles.image); % For complex data
%code under construction .. uncomment when complete / testing - csw
handles.Line = [];
handles.TaskedArea = [];
% Measure Distance Button
handles.measureButton=javaObjectEDT('javax.swing.JToggleButton');
handles.measureButton.setToolTipText('Measure Distance');
handles.measureButton.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
   fullfile(toolboxdir('matlab'),['icons' filesep 'tool_line.png'])));
handles.measureButton.setFocusable(false); % Make consistent with MATLAB toolbars
imageHandle = get(handles.mitm_hand.AxesHandle, 'Children');
%toolbarHandle = get(handles.mitm_hand.AxesHandle, 'Parent');
set(handle(handles.measureButton,'callbackproperties'),'ActionPerformedCallback',@(obj, eventdata) toggleLine(imageHandle));
set(imageHandle, 'ButtonDownFcn',@(obj, eventdata) clickImage(obj));
%Add buttons to toolbar
main_toolbar = handles.mitm_hand.main_toolbar;
main_toolbar.add(handles.measureButton);
main_toolbar.repaint;
main_toolbar.revalidate;
%--end under construction code
handles.mitm_hand.PreChangeViewFcn = @() saveAoiNativeCoords(hObject);
handles.mitm_hand.PostChangeViewFcn = @() restoreAoiLocalCoords(hObject);
handles.phd.axes_handle = []; % For phase history

set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');
set(handles.metaicon,'Color',[51/255 102/255 153/255]);

%disable appropriate controls
set(handles.nx,'enable','off');
set(handles.ny,'enable','off');
set(handles.SelectAOI,'enable','off');
set(handles.SegmentMap,'enable','off');
set(handles.Jump,'enable','off');
set(handles.JumpZoom,'enable','off');
set(handles.Contrast,'enable','off');

axes(handles.logo);
imshow('TaserLogo.png');
set(handles.logo,'XTick',[]);
set(handles.logo,'YTick',[]);
set(handles.logo,'Box','off');

%set background of textboxes/buttons to figure background
Color = get(handles.figure1,'Color');
set(handles.SelectAOIPanel,'BackgroundColor',Color);
set(handles.ForcePow2Check,'BackgroundColor',Color);
set(handles.nxlabel,'BackgroundColor',Color);
set(handles.nylabel,'BackgroundColor',Color);
set(handles.JumpPanel,'BackgroundColor',Color);
set(handles.RunAlgPanel,'BackgroundColor',Color);
set(handles.LatLabel,'BackgroundColor',Color);
set(handles.LonLabel,'BackgroundColor',Color);
set(handles.image,'BackgroundColor','k');
set(handles.txtMeasure,'BackgroundColor',Color);

handles.Color = Color;

handles.meta = [];
handles.PHDData = 0;
handles.AOI_imrect = [];
handles.AOI_native_coords = [];

handles.Prefs = TaserReadPreferences();
if isempty(handles.Prefs)
    handles.Prefs.ImageOverview = 1;
    handles.Prefs.SICDLocation = 'Colocated';    
    handles.Prefs.MaxResolution = 7;
    handles.Prefs.UseOverview = 1;
    handles.Prefs.NumPulses = 1000;
    handles.Prefs.ProcessingPulses = 100;
    handles.Prefs.dBCheck = 1;
    handles.Prefs.ComplexDefaultCheck = 1;
    handles.Prefs.PHDDefaultCheck = 1;
    handles.Prefs.ComplexTemplate = '';
    handles.Prefs.PHDTemplate = '';
end

%set PHD and Complex meta structure settings (either default or read the
%specified files)
if handles.Prefs.ComplexDefaultCheck == 0
   handles.ComplexIconSettings = LoadTemplateFile(handles.Prefs.ComplexTemplate);   
end
if handles.Prefs.ComplexDefaultCheck == 1
    handles.ComplexIconSettings.Line1 = 'IID/Polarization';
    handles.ComplexIconSettings.Line2 = 'Time';
    handles.ComplexIconSettings.Line3 = 'GEO/CC';
    handles.ComplexIconSettings.Line4 = 'IPR';
    handles.ComplexIconSettings.Line5 = 'Graze';
    handles.ComplexIconSettings.Line6 = 'Azimuth';
    handles.ComplexIconSettings.Line7 = 'Layover';
    handles.ComplexIconSettings.Line8 = 'Multipath';
    handles.ComplexIconSettings.Classification = '';
    handles.ComplexIconSettings.Icon = 'true';
    handles.ComplexIconSettings.Direction = 'true';
    handles.ComplexIconSettings.Ascend = 'true';  
end
if handles.Prefs.PHDDefaultCheck == 0
   handles.PHDIconSettings = LoadTemplateFile(handles.Prefs.PHDTemplate);   
end
if handles.Prefs.PHDDefaultCheck == 1
    handles.PHDIconSettings.Line1 = 'IID/Polarization';
    handles.PHDIconSettings.Line2 = 'Time';
    handles.PHDIconSettings.Line3 = 'GEO/CC';
    handles.PHDIconSettings.Line4 = 'IPR';
    handles.PHDIconSettings.Line5 = 'Graze';
    handles.PHDIconSettings.Line6 = 'Azimuth';
    handles.PHDIconSettings.Line7 = 'Layover';
    handles.PHDIconSettings.Line8 = 'Multipath';
    handles.PHDIconSettings.Classification = '';
    handles.PHDIconSettings.Icon = 'true';
    handles.PHDIconSettings.Direction = 'true';
    handles.PHDIconSettings.Ascend = 'true';
end

handles.Algorithms = TaserReadAlgorithms();

%populate list box with algorithms that are selected
UpdateAlgorithmList(handles,handles.Algorithms);

if ~isempty(varargin)
    filename=varargin{1};

    % Update handles structure
    guidata(hObject, handles);

    [path, fname]= fileparts(filename);
    handles.fname = fname;
    handles.NumFiles = 1;
    handles.filename{1} = filename; % Close any old data that was open
    handles = LoadImage(handles);
end


% Update handles structure
guidata(hObject, handles);

function toggleLine(obj)
    handles = guidata(obj);
    if handles.measureButton.isSelected()
        otherButtons = handles.mitm_hand.main_toolbar.getComponents();
        for i = 2:5 % Hardcoded because we know where toggle buttons are        
            otherButtons(i).setSelected(0);
        end
        zoom(handles.figure1,'off');
        pan(handles.figure1,'off');
        datacursormode(handles.figure1,'off');
    else
           set(handles.txtMeasure, 'Visible','off');
           delete(handles.Line);
           handles.Line = [];  
    end
    % Update handles structure
    guidata(obj, handles);
function clickImage(obj)
    handles = guidata(obj);
    %measures the distance.
    if handles.measureButton.isSelected()
        delete(handles.Line);
        guidata(obj, handles);
       [distance bearing handles.Line] = DrawLine(handles, handles.Line);
       set(handles.txtMeasure, 'Visible','on');
       set(handles.txtMeasure, 'String',[num2str(distance) ' m, ' num2str(bearing) ' degrees' ]);
    end
    % Update handles structure
    guidata(obj, handles);

function UpdateAlgorithmList(handles,Algorithms)

count = 0;
for i=1:length(Algorithms)
    if strcmp(Algorithms(i).Selected,'Yes')
        count = count+1;
        AlgList{count} = Algorithms(i).TextName;
    end
end

if (count > 0)
    set(handles.AlgorithmList,'String',AlgList);
else
    set(handles.AlgorithmList,'String',[]);
end

% --- Outputs from this function are returned to the command line.
function varargout = Taser_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
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

%get filename
[fname, path] = uigetfile(sar_file_extensions({'complex','phd'}),'Select Input File',pathstr,'MultiSelect', 'on');   
if isnumeric(fname)
    return;
end

if (iscell(fname))
    NumFiles = length(fname);    
    if (NumFiles > 4) 
        msgbox('Cannot process more than 4 files')
        return;
    end
    for i=1:NumFiles
        filenames(i) = strcat(path,fname(i));
    end           
else
    NumFiles = 1;
    filenames = {strcat(path,fname)};
end

handles.fname = fname;
handles.NumFiles = NumFiles;
handles.filename = filenames; % Close any old data that was open
handles = LoadImage(handles);
setpref('matlab_sar_toolbox','last_used_directory',path); % Store path

% Update handles structure
guidata(hObject, handles);

function handles = LoadImage(handles)

filenames = handles.filename;
NumFiles = handles.NumFiles;
fname = handles.fname;

handles.mitm_hand.close(); % If old data was complex
delete(handles.phd.axes_handle); % If old data was phase history
handles.phd.axes_handle = []; % For phase history
if ~isempty(handles.AOI_imrect)
    delete(handles.AOI_imrect);
    handles.AOI_imrect = [];
    set(handles.SelectAOI,'String','Select AOI');
    guidata(handles.SelectAOI, handles);
end
if ~isempty(handles.TaskedArea)
    delete(handles.TaskedArea);
    handles.TaskedArea = [];
end
guidata(handles.image, handles);

% Open new data
if ~isempty(guess_ph_format(filenames{1})) % Phase history data
    handles.PHDData = 1;
    %determine if we are going to form an image overview or an MPSD
    %first see if we can form an image overview...we must be able to create
    %"cphd" like data if not already passed in, and will only process
    %spotlight mode collects for now
    ph_reader = open_ph_reader(filenames{1});
    phdmeta = ph_reader.get_meta(); 
    if isfield(ph_reader,'read_cphd') && strcmpi(phdmeta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
        FormImage = 1;
    else
        FormImage = 0;
    end
    ph_reader.close();
    if handles.Prefs.ImageOverview && FormImage
        handles.PHDData = 2; %we will allow both complex image, PHD and Image Formation algorithms
        try
            [OverviewName, phdmeta, overviewmeta] = CreateImageOverview(filenames{1},handles.Prefs);
        catch
            errordlg('ERROR: Could Not Create Image Overview.');
            return;
        end
        handles.OverviewName = OverviewName;
        set(handles.image,'BackgroundColor','k');
        set(handles.SelectAOIPanel,'visible','on');
        set(handles.SelectAOI,'enable','on');
        set(handles.SelectAOI,'String','Select AOI');
        set(handles.JumpPanel,'visible','on');
        set(handles.Jump,'enable','on');
        set(handles.JumpZoom,'enable','on');
        set(handles.Contrast,'visible','on');
        set(handles.Contrast,'enable','on');
        handles.mitm_hand.openFile(OverviewName);
        set(handles.mitm_hand.main_toolbar,'Visible',true);    
        meta=handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
        handles.meta = phdmeta;
        handles.overviewmeta = overviewmeta;
        set(handles.SegmentMap,'visible','off');
        if exist('MetaIcon_Custom','file')%may not have this function
            MetaIcon_Custom(meta,handles.ComplexIconSettings,handles.metaicon);   
        else            
            MetaIcon_Complex(meta,'handle',handles.metaicon);            
        end
        
        %we can just add all the algorithms to the list that are selected
        count = 0;
        for i=1:length(handles.Algorithms)
            if strcmp(handles.Algorithms(i).Selected,'Yes') && ...
               ~strcmp(handles.Algorithms(i).Polarimetric,'Yes')
                count = count+1;
                Algorithms{count} = handles.Algorithms(i).TextName; %#ok<AGROW>  
            end
        end
        set(handles.AlgorithmList,'String',Algorithms);
        
        %draw polygon to indicate tasked area
        if exist('DrawTaskedArea.m') > 0
            handles = DrawTaskedArea(handles,phdmeta,overviewmeta);
        end
        
    else    
        set(handles.image,'BackgroundColor',handles.Color);
        set(handles.SelectAOIPanel,'visible','off');
        set(handles.JumpPanel,'visible','off');
        set(handles.SegmentMap,'visible','off');
        set(handles.Contrast,'visible','off');
        set(handles.mitm_hand.main_toolbar,'Visible',false);
        set(handles.mitm_hand.movie_toolbar,'Visible',false);    
        handles = DisplayPHDOverview(handles);
        %Only display PHD (and both) algorithms in Algorithm list    
        count = 0;
        for i=1:length(handles.Algorithms)
            if ~strcmp(handles.Algorithms(i).DataType,'Complex') &&...
               strcmp(handles.Algorithms(i).Selected,'Yes') 
                %add to display list
                count = count+1;
                Algorithms{count} = handles.Algorithms(i).TextName; %#ok<AGROW>
            end
        end
        set(handles.AlgorithmList,'String',Algorithms);    
    end
else % Complex data
    handles.PHDData = 0;
    set(handles.image,'BackgroundColor','k');
    set(handles.SelectAOIPanel,'visible','on');
    set(handles.SelectAOI,'enable','on');
    set(handles.SelectAOI,'String','Select AOI');
    set(handles.JumpPanel,'visible','on');
    set(handles.Jump,'enable','on');
    set(handles.JumpZoom,'enable','on');
    set(handles.Contrast,'visible','on');
    set(handles.Contrast,'enable','on');
    handles.mitm_hand.openFile(filenames);
    set(handles.mitm_hand.main_toolbar,'Visible',true);    
    meta=handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
    handles.meta = meta;
    if (length(handles.mitm_hand.Metadata)>1)
        set(handles.SegmentMap,'visible','on');
        set(handles.SegmentMap,'enable','on');
    else
        set(handles.SegmentMap,'visible','off');
    end
    
    if exist('MetaIcon_Custom','file')%may not have this function
        try 
            MetaIcon_Custom(meta,handles.ComplexIconSettings,handles.metaicon);   
        catch
        end
    else
        try %may not have complete metadata
            MetaIcon_Complex(meta,'handle',handles.metaicon);
        catch
        end
    end
    
    %Only display Complex algorithms in Algorithm list
    %olny display multi-pol algorithms for pol data
    PolData = (numel(handles.mitm_hand.FrameFileIndices{handles.mitm_hand.Frame})>1 ||...
            numel(handles.mitm_hand.FrameSubFileIndices{handles.mitm_hand.Frame})>1);
    count = 0;
    Algorithms = [];
    for i=1:length(handles.Algorithms)
        if ~strcmp(handles.Algorithms(i).DataType,'Phase History') &&...
           strcmp(handles.Algorithms(i).Selected,'Yes') 
            if strcmp(handles.Algorithms(i).Polarimetric,'No') || ...
               strcmp(handles.Algorithms(i).Polarimetric,'All') || PolData 
                %add to display list
                count = count+1;
                Algorithms{count} = handles.Algorithms(i).TextName; %#ok<AGROW>
            end
        end
    end
    set(handles.AlgorithmList,'String',Algorithms);
end

%update Figure Title to include file name(s)
if (NumFiles > 1)
   titlename = fname{1};
else
   titlename = fname;
end
set(handles.figure1,'Name',['Taser, Filename: ',titlename]);

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on button press in SelectAOI.
function SelectAOI_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.AOI_imrect)
    delete(handles.AOI_imrect);
    handles.AOI_imrect = [];
end

if strcmp(get(handles.SelectAOI,'String'),'Clear AOI')
    set(handles.SelectAOI,'String','Select AOI');
else
    %set up imrect in 
    r = getpixelposition(handles.mitm_hand.AxesHandle);
    pos = [r(3)/4 r(4)/4 r(3)/2-1 r(4)/2-1];      
    handles.AOI_imrect = imrect(handles.mitm_hand.AxesHandle,pos);
    zoom(handles.figure1,'off');
    pan(handles.figure1,'off');

    if get(handles.ForcePow2Check,'Value')
        setPositionConstraintFcn(handles.AOI_imrect,@ForcePow2);
    else
        setPositionConstraintFcn(handles.AOI_imrect,@(x) mitm_constraint_fcn(x,handles));
    end

    api_handle = iptgetapi(handles.AOI_imrect);
    api_handle.addNewPositionCallback(@newpos);
    
    newpos(pos);

    set(handles.SelectAOI,'String','Clear AOI');
end

guidata(hObject, handles);

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


function newpos(pos)

%update nx/ny with new values
handles = guidata(gca);
set(handles.nx,'String',round(pos(3)*handles.mitm_hand.Zoom));
set(handles.ny,'String',round(pos(4)*handles.mitm_hand.Zoom));


% --- Executes on button press in ForcePow2Check.
function ForcePow2Check_Callback(hObject, eventdata, handles)
% hObject    handle to ForcePow2Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ForcePow2Check
%if imrect is set, then change the constraint method
if ~isempty(handles.AOI_imrect)
    if get(hObject,'Value')
        setPositionConstraintFcn(handles.AOI_imrect,@ForcePow2);
    else
        setPositionConstraintFcn(handles.AOI_imrect,@(x) mitm_constraint_fcn(x,handles));
    end
end

function pos = ForcePow2(newpos)

handles = guidata(gca);

newsize = newpos(3:4).*handles.mitm_hand.Zoom;
pow2size = (2.^nextpow2(newsize))./handles.mitm_hand.Zoom;
pos = mitm_constraint_fcn([newpos(1:2) pow2size], handles);

function axes_constrained = mitm_constraint_fcn(axes_pos, handles)
    image_pos = round([handles.mitm_hand.axescoords2native(axes_pos(1:2))...
        axes_pos(3:4)*handles.mitm_hand.Zoom]);
    fcn = makeConstrainToRectFcn('imrect',...
        [1 handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageData.NumCols],...
        [1 handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageData.NumRows]);
    image_constrained = double(fcn(image_pos));
    axes_constrained = round([handles.mitm_hand.native2axescoords(image_constrained(1:2))...
        axes_pos(3:4)]);


% --- Executes on button press in MetaViewer.
function MetaViewer_Callback(hObject, eventdata, handles)
% hObject    handle to MetaViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.meta)
    MetaViewer;
else
    MetaViewer(handles.meta);
end


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Message = sprintf('TASER: Version 1.0\nMay 25th, 2011\nTim Cox, Naval Research Lab');
%msgbox(Message,'About TASER','help');
taser_about();

% --------------------------------------------------------------------
function Prefrences_Callback(hObject, eventdata, handles)
% hObject    handle to Prefrences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Prefs = TaserPreferences(handles.Prefs);

%update Icon settings if a file is specified
if ~handles.Prefs.ComplexDefaultCheck
    try
        handles.ComplexIconSettings = LoadTemplateFile(handles.Prefs.ComplexTemplate);   
    catch
        msgbox('Invalid Complex Icon Template File!!');
    end
else
    %reset back to default
    handles.ComplexIconSettings.Line1 = 'IID/Polarization';
    handles.ComplexIconSettings.Line2 = 'Time';
    handles.ComplexIconSettings.Line3 = 'GEO/CC';
    handles.ComplexIconSettings.Line4 = 'IPR';
    handles.ComplexIconSettings.Line5 = 'Graze';
    handles.ComplexIconSettings.Line6 = 'Azimuth';
    handles.ComplexIconSettings.Line7 = 'Layover';
    handles.ComplexIconSettings.Line8 = 'Multipath';
    handles.ComplexIconSettings.Classification = '';
    handles.ComplexIconSettings.Icon = 'true';
    handles.ComplexIconSettings.Direction = 'true';
    handles.ComplexIconSettings.Ascend = 'true';  
end

if ~handles.Prefs.PHDDefaultCheck
    try
        handles.PHDIconSettings = LoadTemplateFile(handles.Prefs.PHDTemplate);   
    catch
        msgbox('Invalid PHD Icon Template File!!');
    end
else
    %reset back to default  
    handles.PHDIconSettings.Line1 = 'IID/Polarization';
    handles.PHDIconSettings.Line2 = 'Time';
    handles.PHDIconSettings.Line3 = 'GEO/CC';
    handles.PHDIconSettings.Line4 = 'IPR';
    handles.PHDIconSettings.Line5 = 'Graze';
    handles.PHDIconSettings.Line6 = 'Azimuth';
    handles.PHDIconSettings.Line7 = 'Layover';
    handles.PHDIconSettings.Line8 = 'Multipath';
    handles.PHDIconSettings.Classification = '';
    handles.PHDIconSettings.Icon = 'true';
    handles.PHDIconSettings.Direction = 'true';
    handles.PHDIconSettings.Ascend = 'true';
end
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Options_Callback(hObject, eventdata, handles)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Algorithms_Callback(hObject, eventdata, handles)
% hObject    handle to Algorithms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Algorithms = AlgorithmSelection(handles.Algorithms);

UpdateAlgorithmList(handles,handles.Algorithms)

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in AlgorithmList.
function AlgorithmList_Callback(hObject, eventdata, handles)
% hObject    handle to AlgorithmList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AlgorithmList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AlgorithmList


% --- Executes during object creation, after setting all properties.
function AlgorithmList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlgorithmList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get Algorithm name from Listbox
ListVal = get(handles.AlgorithmList,'Value');
ListStrings = get(handles.AlgorithmList,'String');
Algorithm = ListStrings{ListVal};

%loop through Algorithm Structure to determine Run Name.  Also, provide
%some checking for algorithm type (Data, polarimetric etc.)
for i=1:length(handles.Algorithms)
    if strcmp(Algorithm,handles.Algorithms(i).TextName)
        CallName = handles.Algorithms(i).CallName;
        DataType = handles.Algorithms(i).DataType;
        Polarimetric = handles.Algorithms(i).Polarimetric;
        break;
    end
end

%Make sure algorithm runs the type of dta we have (PHD/Complex)
if (handles.PHDData == 1 && strcmpi(DataType,'Complex'))
    msgbox('Algorithm requires Complex Data!');
    return;
end

if (handles.PHDData == 0 && strcmpi(DataType,'Phase History'))
    msgbox('Algorithm requires Phase History Data!');
    return;
end

%if no file has been loaded then just launch tool
if isempty(handles.meta)
    eval(CallName);
    return;
end

%get AOI
if ~isempty(handles.AOI_imrect)
    axes_pos = getPosition(handles.AOI_imrect);
    aoi = round([handles.mitm_hand.axescoords2native(axes_pos(1:2))...
        axes_pos(3:4)*handles.mitm_hand.Zoom]);
    im_data = handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageData; % For ease of notation
    aoi(3:4) = min(double([im_data.NumCols im_data.NumRows])-aoi(1:2),aoi(3:4));
else % Otherwise use entire visible area
    aoi = handles.mitm_hand.Position;
end

%make sure this is the right type of data
%todo: check data type (complex, phd) once phd data is loaded into tool

if handles.PHDData == 1
    filenames = handles.filename;
    feval(CallName, 'filename', filenames{1});
    return;   
end

if handles.PHDData == 2
    %call based on algorithm type
    if strcmpi(DataType,'Complex') || strcmpi(DataType,'All')
        feval(CallName, 'filename', handles.OverviewName, ...
            'aoi', aoi, 'segment', handles.mitm_hand.Frame);
    elseif strcmpi(DataType,'Processing')
        %get center lat/lon, range/az size (slant)
        CenterX = aoi(1) + aoi(3)/2;
        CenterY = aoi(2) + aoi(4)/2;
        pos = point_slant_to_ground([CenterY;CenterX], handles.overviewmeta); 
        RangeSize = aoi(3)*handles.overviewmeta.Grid.Row.SS;
        AzSize = aoi(4)*handles.overviewmeta.Grid.Col.SS;
        filenames = handles.filename;
        feval(CallName, 'filename', filenames{1},'center',pos,'size',[AzSize RangeSize]);
    else
        filenames = handles.filename;
        feval(CallName, 'filename', filenames{1});
    end
    return;
end

filenames = handles.mitm_hand.Filename(handles.mitm_hand.FrameFileIndices{handles.mitm_hand.Frame});
segments = handles.mitm_hand.FrameSubFileIndices{handles.mitm_hand.Frame};
if numel(filenames)>1||numel(segments)>1 % Data is multi-pol
    DataPol = 1;
else
    DataPol = 0;
end

if (DataPol==0 &&  strcmp(Polarimetric,'Yes'))
    choice = questdlg('Algorithm requires multi-pol data.  Do you want to launch?',...
                      'Warning','No','Yes','No');
    if strcmp(choice,'No')
        return;
    else
        %run without any arguments
        eval(CallName);
        return;
    end
end

if (DataPol==1 && strcmp(Polarimetric,'No'))
    Polarizations = handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageFormation.TxRcvPolarizationProc;
    if ~iscell(Polarizations)
        Polarizations = {Polarizations};
    end
    %allow user to select data channel to run (VV,HH,HV or VH) or to cancel
    %all together
    if (length(Polarizations) == 2) %#ok<*ISMT>
        choice = questdlg('Algorithm only processes one data channel. Process single channel?', ...
            'Warning',Polarizations{1},Polarizations{2},'No','No');
    elseif (length(Polarizations) == 3)
        choice = questdlg('Algorithm only processes one data channel. Process single channel?', ...
            'Warning',Polarizations{1},Polarizations{2},...
            Polarizations{3},'No','No');
    elseif (length(Polarizations) == 4)
        choice = questdlg('Algorithm only processes one data channel. Process single channel?', ...
            'Warning',Polarizations{1},Polarizations{2},...
            Polarizations{3},Polarizations{4},'No','No');
    end
    
    %now call algorithm with specified filename
    if (strcmp(choice,Polarizations{1}))
        feval(CallName, 'filename', handles.filename{1}, ...
            'aoi', aoi, 'segment', segments(1));
    elseif (strcmp(choice,Polarizations{2}))
        feval(CallName, 'filename', handles.filename{2}, ...
            'aoi', aoi, 'segment', segments(2));
    elseif ((length(Polarizations) > 2) && strcmp(choice,Polarizations{3}))
        feval(CallName, 'filename', handles.filename{3}, ...
            'aoi', aoi, 'segment', segments(3));
    elseif ((length(Polarizations) > 3) && strcmp(choice,Polarizations{4}))
        feval(CallName, 'filename', handles.filename{4}, ...
            'aoi', aoi, 'segment', segments(4));
    else
        return;
    end
    return;
end

if size(filenames,2) == 1
    feval(CallName, 'filename', filenames{1}, ... % Run tool
    'aoi', aoi, 'segment', segments);
else
    feval(CallName, 'filename', filenames, ... % Run tool
    'aoi', aoi, 'segment', segments);
end


function handles = DisplayPHDOverview(handles)

%determine if PHDExlorer mat file exists.  If so then use the MPSD stored
%in the file.  It is better, and we don't have to compute anything
[pathstr, name, ext] = fileparts(handles.filename{1});
matname = [pathstr filesep name '_Channel_' num2str(1) '.mat'];
listing = dir(matname);

if ~isempty(listing)
    %load mat file and set freq and mpsddata
    load(matname);
    center_vector = round(length(metaicon_args{2}.TxTime)/2);
    SCP = metaicon_args{2}.SRPPos(center_vector,:).';
    lla = ecf_to_geodetic(SCP);
    IID = meta.CollectionInfo.CoreName;
else
    %get phd metadata
    ph_reader_object = open_ph_reader(handles.filename{1});
    meta = ph_reader_object.get_meta();
    NumPulses = meta.Data.ArraySize.NumVectors;
    [ignore, vbmeta] = ph_reader_object.read_raw([1 round(meta.Data.ArraySize(1).NumVectors/2) meta.Data.ArraySize(1).NumVectors],[]);
    IID = meta.CollectionInfo.CoreName;
    center_vector = round(length(vbmeta.TxTime)/2);
    SCP = vbmeta.SRPPos(center_vector,:).';
    lla = ecf_to_geodetic(SCP);
    metaicon_args = {meta, vbmeta};
end

NumMPSDPulses = handles.Prefs.NumPulses;

if ~isempty(listing)
    %load mat file and set freq and mpsddata
    load(matname);
    %get array of freqs
    sampling_rate =  meta.RadarCollection.Waveform.WFParameters.ADCSampleRate;
    chirp = strcmpi('CHIRP',meta.RadarCollection.Waveform.WFParameters.RcvDemodType);
    if chirp
        fc = meta.RadarCollection.Waveform.WFParameters.RcvFreqStart;
        freq = linspace(fc-sampling_rate/2,fc+sampling_rate/2,meta.Data.ArraySize(1).NumSamples);
    else
        deramp_rate = meta.RadarCollection.Waveform.WFParameters.RcvFMRate;
        freq = meta.RadarCollection.Waveform.WFParameters.RcvFreqStart + ...
            ((0:(double(meta.Data.ArraySize(1).NumSamples)-1))' * deramp_rate / sampling_rate);
    end
    %convert pulses to db
    sum_pulses(sum_pulses <= 0) = 1;
    sum_pulses = 10.*log10(sum_pulses);
    mpsddata = sum_pulses - min(sum_pulses); % Scale such that min is at 0 dB
else
    %compute from phd file
    [freq mpsddata] = MPSD(ph_reader_object,0,1,NumPulses,floor(NumPulses/NumMPSDPulses),handles.Prefs.ProcessingPulses);
end

%now display a MPSD plot and map with the image location marked
handles.phd.axes_handle(1) = axes('Parent',handles.image,'Position',[.075 .075 .85 .4]);
plot(freq./10^9,mpsddata);
xlabel('Frequency (GHz)');
ylabel('Power (dB)');
xlim([min(freq./10^9) max(freq./10^9)]); 
title(['Mean Power Spectral Density for IID: ' IID]);
handles.phd.axes_handle(2) = axes('Parent',handles.image,'Position',[.075 .55 .85 .4]);
%plot coastal map
[Lats Lons] = GetCoastalPoints(1);
hold on;
plot(Lons,Lats,'-k');
ylim([-90 90]);
xlim([-180 180]);
grid on; box on;
%plot target point
plot(lla(2),lla(1),'MarkerSize',10,'Marker','o','MarkerEdgeColor','Red','MarkerFaceColor','Red');
hold off;

%display metaicon
MetaIcon_Custom(meta,handles.PHDIconSettings,handles.metaicon); 

if isempty(listing)
    ph_reader_object.close();
end
handles.meta = meta;


% --- Executes on button press in SegmentMap.
function SegmentMap_Callback(hObject, eventdata, handles)
% hObject    handle to SegmentMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if iscell(handles.filename)
    fname = handles.filename{1};
else
    fname = handles.filename;
end

reader_obj = open_reader(fname);

if ~iscell(reader_obj)
    reader_obj{1} = reader_obj;
end

figure;
hold on;

for i=1:length(reader_obj)
    meta = reader_obj{i}.get_meta();
    ICP = meta.GeoData.ImageCorners.ICP;
    Lats = [ICP.FRFC.Lat ICP.FRLC.Lat ICP.LRLC.Lat ICP.LRFC.Lat ICP.FRFC.Lat];
    Lons = [ICP.FRFC.Lon ICP.FRLC.Lon ICP.LRLC.Lon ICP.LRFC.Lon ICP.FRFC.Lon];
    text(mean(Lons(1:4)),mean(Lats(1:4)),num2str(i));
    plot(Lons,Lats);    
end

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

IID = meta.CollectionInfo.CoreName(1:min(16,end));
a = sprintf('Segment Map for Collect: %s',IID);
title(a);

hold off;


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function Latitude_Callback(hObject, eventdata, handles)
% hObject    handle to Latitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Latitude as text
%        str2double(get(hObject,'String')) returns contents of Latitude as a double


% --- Executes during object creation, after setting all properties.
function Latitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Latitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Longitude_Callback(hObject, eventdata, handles)
% hObject    handle to Longitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Longitude as text
%        str2double(get(hObject,'String')) returns contents of Longitude as a double


% --- Executes during object creation, after setting all properties.
function Longitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Longitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Jump.
function Jump_Callback(hObject, eventdata, handles)
% hObject    handle to Jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get lat/lon and convert to decimal degrees
LatStr = get(handles.Latitude,'String');
LonStr = get(handles.Longitude,'String');

%determine if both strings are in Lat field...we will only split this if
%there is a N or S
Index1 = strfind(LatStr,'N');
Index2 = strfind(LatStr,'S');
if ~isempty(Index1) && isnan(latlonnum(LonStr))
    LonStr = LatStr(Index1+1:end);
    LatStr = LatStr(1:Index1);
end
if ~isempty(Index2) && isnan(latlonnum(LonStr))
    LonStr = LatStr(Index2+1:end);
    LatStr = LatStr(1:Index2);
end

Lat = latlonnum(LatStr);
Lon = latlonnum(LonStr);
try
    Alt = handles.meta.GeoData.SCP.LLH.HAE;
catch
    Alt = 0;
end
lla = [Lat;Lon;Alt];

Replaced = 0;

TotSegments = numel(handles.mitm_hand.Metadata);
if TotSegments > 1
    if iscell(handles.filename)
        fname = handles.filename{1};
    else
        fname = handles.filename;
    end
    reader_obj = open_reader(fname);
    count = 0;
    for i=1:TotSegments
        meta = reader_obj{i}.get_meta();
        row_col_pos = point_ground_to_slant(lla, meta, 'projectToDEM', false).';
        XPos = row_col_pos(2); YPos = row_col_pos(1);
        nx = meta.ImageData.NumCols;
        ny = meta.ImageData.NumRows;
        if (XPos > 1) && (XPos < nx) && (YPos > 1) && (YPos < ny)
            count = count + 1;
            %store segment
            Segments(count) = i;
            %compute minimum distance to edge.  If the point is in multiple
            %segments then we'll use the segment in which the point is the
            %farthest away from any edge
            d1 = XPos;
            d2 = YPos;
            d3 = nx - XPos;
            d4 = ny - YPos;
            MinDist(count) = min([d1 d2 d3 d4]);
        end
    end
    %set segment and meta and XPos and YPos
    [val index] = max(MinDist);
    handles.mitm_handle.Frame = Segments(index);
else
    %check to see if user passed in X/Y pos.  We'll assume it's an X/Y
    %position if the lat/lon is not in the image, the X and Y values are
    %integers, and the position X/Y position is in the image
    meta = handles.mitm_hand.Metadata{1};
    nx = meta.ImageData.NumCols;
    ny = meta.ImageData.NumRows;
    pos = point_ground_to_slant(lla,meta);
    X = pos(2); Y = pos(1);
    if X > 1 && X <= nx && Y > 1 && Y <=ny
        InImage = 1;
    else
        InImage = 0;
    end
    xpos = str2double(LatStr);
    ypos = str2double(LonStr);
    if round(xpos) == xpos && round(ypos) == ypos && xpos > 1 && ...
       xpos <= nx && ypos > 1 && ypos <= ny && ~InImage  
        %compute lla from pixel position
        lla = point_slant_to_ground([ypos;xpos;Alt], meta);
        Replaced = 1;
    end
end
    
if ~Replaced
    if isnan(Lat)
        msgbox('Invalid Latitude!!');
        return;
    end

    if isnan(Lon)
        msgbox('Invalid Longitude!!');
        return;
    end
end

if ~handles.mitm_hand.geojump(lla)
    msgbox('Specified point is not in image!');
end

% Update handles structure
guidata(hObject, handles);

function saveAoiNativeCoords(src)

handles = guidata(src);

if ~isempty(handles.AOI_imrect)
    local_coords = handles.AOI_imrect.getPosition();
    native_coords = handles.mitm_hand.axescoords2native(local_coords(1:2));
    native_coords(3:4) = local_coords(3:4)*handles.mitm_hand.Zoom;
    handles.AOI_native_coords = native_coords;
else
    handles.AOI_native_coords = [];
end
if ~isempty(handles.Line)
    local_coords = handles.Line.getPosition();
    native_coords(1,:) = handles.mitm_hand.axescoords2native(local_coords(1,:));
    native_coords(2,:) = handles.mitm_hand.axescoords2native(local_coords(2,:));
    handles.Line_native_coords = native_coords;
else
    handles.Line_native_coords = [];
end
if ~isempty(handles.TaskedArea)
    XData = get(handles.TaskedArea,'XData');
    YData = get(handles.TaskedArea,'YData');
    local_coords = [XData' YData'];
    native_coords(1,:) = handles.mitm_hand.axescoords2native(local_coords(1,:));
    native_coords(2,:) = handles.mitm_hand.axescoords2native(local_coords(2,:));
    native_coords(3,:) = handles.mitm_hand.axescoords2native(local_coords(3,:));
    native_coords(4,:) = handles.mitm_hand.axescoords2native(local_coords(4,:));
    native_coords(5,:) = handles.mitm_hand.axescoords2native(local_coords(5,:));
    handles.TaskedArea_native_coords = native_coords;
else
    handles.TaskedArea_native_coords = [];
end


guidata(src,handles);
    
function restoreAoiLocalCoords(src)

oldzoom = get(zoom(src),'Enable');
oldpan = get(pan(src),'Enable');

handles = guidata(src);
if ~isempty(handles.AOI_native_coords)
    delete(handles.AOI_imrect);
    local_coords = handles.mitm_hand.native2axescoords(...
        handles.AOI_native_coords(1:2));
    local_coords(3:4) = handles.AOI_native_coords(3:4)/handles.mitm_hand.Zoom;
    handles.AOI_imrect = imrect(handles.mitm_hand.AxesHandle, local_coords);
    handles.AOI_native_coords = [];
    % Restore callbacks to imrect
    api_handle = iptgetapi(handles.AOI_imrect);
    api_handle.addNewPositionCallback(@newpos);
    ForcePow2Check_Callback(handles.ForcePow2Check, [], handles);
end
if ~isempty(handles.Line_native_coords)
    delete(handles.Line);
    
    local_coords(1,:) = handles.mitm_hand.native2axescoords(...
        handles.Line_native_coords(1,:)); 
    local_coords(2,:) = handles.mitm_hand.native2axescoords(...
        handles.Line_native_coords(2,:)); 
    handles.Line = imline(handles.mitm_hand.AxesHandle, local_coords);
    handles.Line_native_coords = [];
    setColor(handles.Line,'Cyan');
    kids = get(handles.Line,'Children');
    set(kids,'LineWidth',2);
    set(kids,'MarkerSize',1);
    %add callback for user moving line
    api_handle = iptgetapi(handles.Line);
    api_handle.addNewPositionCallback(@newlinepos);
end
if ~isempty(handles.TaskedArea_native_coords)
    delete(handles.TaskedArea);
    local_coords(1,:) = handles.mitm_hand.native2axescoords(...
        handles.TaskedArea_native_coords(1,:)); 
    local_coords(2,:) = handles.mitm_hand.native2axescoords(...
        handles.TaskedArea_native_coords(2,:)); 
    local_coords(3,:) = handles.mitm_hand.native2axescoords(...
        handles.TaskedArea_native_coords(3,:)); 
    local_coords(4,:) = handles.mitm_hand.native2axescoords(...
        handles.TaskedArea_native_coords(4,:)); 
    local_coords(5,:) = handles.mitm_hand.native2axescoords(...
        handles.TaskedArea_native_coords(5,:)); 
    
    %set(gcf,'CurrentAxes',handles.mitm_hand.AxesHandle);    
    hold(handles.mitm_hand.AxesHandle,'on');
    handles.TaskedArea = plot(handles.mitm_hand.AxesHandle,local_coords(:,1),local_coords(:,2),'-r','LineWidth',2);    
    hold(handles.mitm_hand.AxesHandle,'off');
end

guidata(src,handles);

function newlinepos(pos)

handles = guidata(gcbo);
distance = ComputeDistance(handles.mitm_hand, pos);
set(handles.txtMeasure, 'String', [num2str(distance) ' m']);

% Restore old state which was changed when we redrew the shapes
try % Doesn't work when initiated from zoom callback.
    zoom(src,oldzoom);
    pan(src,oldpan);
end


% --- Executes on button press in Contrast.
function Contrast_Callback(hObject, eventdata, handles)
% hObject    handle to Contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.mitm_hand.Remap = 'linearremap';
imcontrast(handles.mitm_hand.AxesHandle);


% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


% --- Executes on button press in JumpZoom.
function JumpZoom_Callback(hObject, eventdata, handles)
% hObject    handle to JumpZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get lat/lon and convert to decimal degrees
LatStr = get(handles.Latitude,'String');
LonStr = get(handles.Longitude,'String');

%determine if both strings are in Lat field...we will only split this if
%there is a N or S
Index1 = strfind(LatStr,'N');
Index2 = strfind(LatStr,'S');
if ~isempty(Index1) && isnan(latlonnum(LonStr))
    LonStr = LatStr(Index1+1:end);
    LatStr = LatStr(1:Index1);
end
if ~isempty(Index2) && isnan(latlonnum(LonStr))
    LonStr = LatStr(Index2+1:end);
    LatStr = LatStr(1:Index2);
end

Lat = latlonnum(LatStr);
Lon = latlonnum(LonStr);
try
    Alt = handles.meta.GeoData.SCP.LLH.HAE;
catch
    Alt = 0;
end
lla = [Lat;Lon;Alt];

Replaced = 0;

TotSegments = numel(handles.mitm_hand.Metadata);
if TotSegments > 1
    if iscell(handles.filename)
        fname = handles.filename{1};
    else
        fname = handles.filename;
    end
    reader_obj = open_reader(fname);
    count = 0;
    for i=1:TotSegments
        meta = reader_obj{i}.get_meta();
        row_col_pos = point_ground_to_slant(lla, meta, 'projectToDEM', false).';
        XPos = row_col_pos(2); YPos = row_col_pos(1);
        nx = meta.ImageData.NumCols;
        ny = meta.ImageData.NumRows;
        if (XPos > 1) && (XPos < nx) && (YPos > 1) && (YPos < ny)
            count = count + 1;
            %store segment
            Segments(count) = i;
            %compute minimum distance to edge.  If the point is in multiple
            %segments then we'll use the segment in which the point is the
            %farthest away from any edge
            d1 = XPos;
            d2 = YPos;
            d3 = nx - XPos;
            d4 = ny - YPos;
            MinDist(count) = min([d1 d2 d3 d4]);
        end
    end
    %set segment and meta and XPos and YPos
    [val index] = max(MinDist);
    handles.mitm_handle.Frame = Segments(index);
else
    %check to see if user passed in X/Y pos.  We'll assume it's an X/Y
    %position if the lat/lon is not in the image, the X and Y values are
    %integers, and the position X/Y position is in the image
    meta = handles.mitm_hand.Metadata{1};
    nx = meta.ImageData.NumCols;
    ny = meta.ImageData.NumRows;
    pos = point_ground_to_slant(lla,meta);
    X = pos(2); Y = pos(1);
    if X > 1 && X <= nx && Y > 1 && Y <=ny
        InImage = 1;
    else
        InImage = 0;
    end
    xpos = str2double(LatStr);
    ypos = str2double(LonStr);
    if round(xpos) == xpos && round(ypos) == ypos && xpos > 1 && ...
       xpos <= nx && ypos > 1 && ypos <= ny && ~InImage  
        %compute lla from pixel position
        lla = point_slant_to_ground([ypos;xpos;Alt], meta);
        Replaced = 1;
    end
end
    
if ~Replaced
    if isnan(Lat)
        msgbox('Invalid Latitude!!');
        return;
    end

    if isnan(Lon)
        msgbox('Invalid Longitude!!');
        return;
    end
end

if InImage
    handles.mitm_hand.setView('CenterPos',pos,'Zoom',1);
end

if ~handles.mitm_hand.geojump(lla)    
    msgbox('Specified point is not in image!');
end

% Update handles structure
guidata(hObject, handles);
