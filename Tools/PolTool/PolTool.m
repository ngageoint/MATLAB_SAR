function varargout = PolTool(varargin)
%PolTool Interactive Polarization Combination Tool
%
% Applies weighted complex combination of linear polarized data.  For
% Quad-Pol this allows us to synthesize any transmit and recieve
% polarimetric configuration.  Dual Pol allows us to only synthesize one 
% side (typically Rx). 
% Polarization is defined by angle and ellipticity.  Angle is defined from
% 0-360 (where 0 is H).  Ellipticity is defined from -45 to 45.  Positive
% ellipticity indicated CW circular rotation.
% Interactive control is performed with the mouse. Motion with the left
% button down modifies the Tx polarization and right button down modifies
% the Rx polarization.  Horizontal motion changes the angle (Right is
% positive), vertical motion changes the ellipticity (Up is positive). The
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
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

gui_Singleton = 0;
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

% Setup image display
handles.mitm_hand = hg_mitm_viewer(handles.image);
% Could use callbacks to update pixel/AOI selections to sync with pan/zoom
% handles.mitm_hand.PreChangeViewFcn = @() saveShapesNativeCoords(hObject);
% handles.mitm_hand.PostChangeViewFcn = @() restoreShapesLocalCoords(hObject);

%format metaicon
set(handles.metaicon,'XTick',[]);
set(handles.metaicon,'YTick',[]);
set(handles.metaicon,'Box','on');
set(handles.metaicon,'Color',[51/255 102/255 153/255]);
setAllowAxesZoom(zoom(handles.figure1),handles.metaicon,false);
setAllowAxesPan(pan(handles.figure1),handles.metaicon,false);

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

%set defaults
handles.TxAng = 0;
handles.RxAng = 0;
handles.TxEllip = 0;
handles.RxEllip = 0;
handles.TxAngleFloat = 0;
handles.RxAngleFloat = 0;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
handles.Select = 'None';
handles.InMotionFcn = 0;
handles.SelectingRegion = 0;
handles.aviflag = 0;
handles.PixelSelect = 0;
handles.PointInImage = 0;
handles.XOffset = 0;
handles.YOffset = 0;
handles.aoi = [];
set(handles.AngleControl,'Value',1);
set(handles.EllipControl,'Value',1);
set(handles.TxAngleInc,'String',1);
set(handles.RxAngleInc,'String',1);
set(handles.TxEllipInc,'String',1);
set(handles.RxEllipInc,'String',1);
set(handles.TxAngleMan,'String',0);
set(handles.RxAngleMan,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
set(handles.FrameRate,'String',5);
set(handles.Grayscale,'Value',1);
set(handles.AngRes,'String',5);

% Update handles structure
guidata(hObject, handles);

set(handles.TxAngle,'String',handles.TxAng);
set(handles.RxAngle,'String',handles.RxAng);
set(handles.TxEllipticity,'String',handles.TxEllip);
set(handles.RxEllipticity,'String',handles.RxEllip);

%add radio-button callbacks for color/grayscale
set(handles.Display,'SelectionChangeFcn',@Colormap_Callback);

%enable/disable controls
set(handles.TxHCheck,'Enable','off');
set(handles.RxHCheck,'Enable','off');
set(handles.TxVCheck,'Enable','off');
set(handles.RxVCheck,'Enable','off');
set(handles.HHButton,'Enable','off');
set(handles.HVButton,'Enable','off');
set(handles.VHButton,'Enable','off');
set(handles.VVButton,'Enable','off');
set(handles.CircCoButton,'Enable','off');
set(handles.CircCrossButton,'Enable','off');
set(handles.HHPlusVVButton,'Enable','off');
set(handles.HHMinusVVButton,'Enable','off');
set(handles.Color,'Enable','off');
set(handles.Grayscale,'Enable','off');
set(handles.SaveImage,'Enable','off');
set(handles.StartMovie,'Enable','off');
set(handles.StopMovie,'Enable','off');
set(handles.RunDecomp,'Enable','off');
set(handles.Generate,'Enable','off');
set(handles.SaveDecomp,'Enable','off');
set(handles.RedPol,'Enable','inactive');
set(handles.GreenPol,'Enable','inactive');
set(handles.BluePol,'Enable','inactive');
set(handles.SetRed,'Enable','off');
set(handles.SetGreen,'Enable','off');
set(handles.SetBlue,'Enable','off');
set(handles.SaveCov,'Enable','off');

set(handles.SelectPixel,'Enable','off');
set(handles.SelectAOI,'Enable','off');
set(handles.PlotResults,'Enable','off');

set(handles.TxManPanel,'Visible','off');
set(handles.RxManPanel,'Visible','off');

%plot Tx and Rx Ellipses
PlotPolEllipse(handles.TxPlot,handles.TxAng,handles.TxEllip);
PlotPolEllipse(handles.RxPlot,handles.RxAng,handles.RxEllip);

pathstr=fileparts(mfilename('fullpath'));

% if isdeployed
%     %Trim off PolTool from path since decomposition folder is at the same
%     %level (not sure why)
%     pathlength = length(pathstr);
%     pathstr = pathstr(1:(pathlength-7));
% end

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

set(handles.DecompList,'String',stringtemp);

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[]);
p.addParamValue('segment',1);
p.parse(varargin{:});

% determine if image name was passed as calling argument

if ~isempty(p.Results.filename)
    LoadImage(p.Results.filename,hObject,handles,p.Results.aoi,p.Results.segment);
else
    % Update handles structure
    guidata(hObject, handles);
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


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isfield(handles,'Select')&&strcmp(handles.Select,'None'))
    return;
end

if (~isfield(handles,'InMotionFcn') || handles.InMotionFcn == 1 || ...
        isempty(handles.mitm_hand.Metadata) || ...
        handles.SelectingRegion == 1 || handles.PixelSelect == 1)
    return;
end

VVal = get(handles.TxVCheck,'Value');
HVal = get(handles.TxHCheck,'Value');
if (strcmp(handles.Select,'Tx') && (VVal == 0 || HVal == 0))
    return;
end

handles.InMotionFcn = 1;
% Update handles structure
guidata(hObject, handles);

AngleControl = get(handles.AngleControl,'Value');
EllipControl = get(handles.EllipControl,'Value');

AngleLock = get(handles.AngleLock,'Value');
EllipLock = get(handles.EllipLock,'Value');

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
set(handles.TxAngleMan,'String',TxAng);
set(handles.TxEllipMan,'String',TxEllip);
set(handles.RxAngleMan,'String',RxAng);
set(handles.RxEllipMan,'String',RxEllip);

handles.pos = pos;

%update Ellipse Plots (only if necessary)
if (strcmp(handles.Select,'Tx') || AngleLock || EllipLock)
    PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);
end
if (strcmp(handles.Select,'Rx') || AngleLock || EllipLock)
    PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);
end

%Update Image
handles.RxAngleFloat = RxAng;
handles.RxEllipticityFloat = RxEllip;
handles.TxAngleFloat = TxAng;    
handles.TxEllipticityFloat = TxEllip;
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
    %Left Mouse Click
    handles.Select = 'Tx';
elseif strcmpi(typ,'alt')
    handles.Select = 'Rx';
else
    if ( get(handles.TxVCheck,'Value') == 0)
        %set to HH
        handles.TxAng = 0;
        handles.TxEllip = 0;
        handles.RxAng = 0;
        handles.RxEllip = 0;
    elseif (get(handles.TxHCheck,'Value') == 0)
        %set to VV
        handles.TxAng = 90;
        handles.TxEllip = 0;
        handles.RxAng = 90;
        handles.RxEllip = 0;
    end    
    %reset all to "default", which varies based on input type
    set(handles.TxAngle,'String',round(handles.TxAng));
    set(handles.TxEllipticity,'String',round(handles.TxEllip));
    set(handles.RxAngle,'String',round(handles.RxAng));
    set(handles.RxEllipticity,'String',round(handles.RxEllip));
    set(handles.TxAngleMan,'String',handles.TxAng);
    set(handles.TxEllipMan,'String',handles.TxEllip);
    set(handles.RxAngleMan,'String',handles.RxAng);
    set(handles.RxEllipMan,'String',handles.RxEllip);    
    handles.TxAngleFloat = handles.TxAng;
    handles.RxAngleFloat = handles.RxAng;
    handles.TxEllipticityFloat = handles.TxEllip;
    handles.RxEllipticityFloat = handles.RxEllip;
    PlotPolEllipse(handles.TxPlot,handles.TxAng,handles.TxEllip);
    PlotPolEllipse(handles.RxPlot,handles.RxAng,handles.RxEllip);
    UpdateImage(handles);
    handles.Select = 'None';
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


function PlotPolEllipse(handle,Ang,Ellip)

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


function Filename_Callback(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Filename as text
%        str2double(get(hObject,'String')) returns contents of Filename as a double


% --- Executes during object creation, after setting all properties.
function Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
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


% --- Executes on button press in HHCheck.
function HHCheck_Callback(hObject, eventdata, handles)
% hObject    handle to HHCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HHCheck


% --- Executes on button press in HVCheck.
function HVCheck_Callback(hObject, eventdata, handles)
% hObject    handle to HVCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HVCheck


% --- Executes on button press in VVCheck.
function VVCheck_Callback(hObject, eventdata, handles)
% hObject    handle to VVCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VVCheck


% --- Executes on button press in VHCheck.
function VHCheck_Callback(hObject, eventdata, handles)
% hObject    handle to VHCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of VHCheck


function LoadImage(iqfiles,hObject,handles,aoi,segment)

% Some file formats allow multiple polarimetric channels to be stored
% within a single file.  Some file formats have polarimetric datasets
% stored across multiple files.

% Input arguments
if ischar(iqfiles) % Input can be single string or cell array of strings
    iqfiles = {iqfiles}; % Just treat both types as cell array
end
set(handles.Filename,'String',iqfiles{1});

% Open file(s)
handles.mitm_hand.DataTransformFcn = [];
handles.mitm_hand.close();
if isempty(aoi) % Default is to show entire image
    handles.mitm_hand.openFile(iqfiles);
else % AOI was passed in
    handles.mitm_hand.openFile(iqfiles, true);
    pixels_available = floor(getpixelposition(handles.mitm_hand.AxesHandle))-1;
    handles.mitm_hand.setView('CenterPos',aoi(1:2)+(aoi(3:4)/2),...
        'Zoom',max(ceil(aoi(3:4)./pixels_available(3:4))),'Frame',segment(1));
end
meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
delete(get(handles.metaicon,'Children')); % Clear old metaicon
try
    MetaIcon_Complex(meta,'handle',handles.metaicon);
end

% Set GUI fields based on data
handles.HH=find(strcmpi('H:H',meta.ImageFormation.TxRcvPolarizationProc));
handles.HV=find(strcmpi('H:V',meta.ImageFormation.TxRcvPolarizationProc));
handles.VH=find(strcmpi('V:H',meta.ImageFormation.TxRcvPolarizationProc));
handles.VV=find(strcmpi('V:V',meta.ImageFormation.TxRcvPolarizationProc));
set(handles.TxHCheck,'Value',~isempty(handles.HH)||~isempty(handles.HV));
set(handles.TxVCheck,'Value',~isempty(handles.VV)||~isempty(handles.VH));
set(handles.RxHCheck,'Value',~isempty(handles.HH)||~isempty(handles.VH));
set(handles.RxVCheck,'Value',~isempty(handles.VV)||~isempty(handles.HV));
set([handles.TxHCheck handles.TxVCheck handles.RxHCheck handles.RxVCheck],'Enable','off');
if ~isempty(handles.HH)
    set(handles.TxHCheck,'Enable','on');
end
if ~isempty(handles.HV)
    set(handles.TxHCheck,'Enable','on');
end
if ~isempty(handles.VH)
    set(handles.TxVCheck,'Enable','on');
end
if ~isempty(handles.VV)
    set(handles.TxVCheck,'Enable','on');
end
if ~isempty(handles.HH)&&~isempty(handles.VV)
    % Complete "tri-pol" datasets
    if isempty(handles.HV)&&~isempty(handles.VH)
        handles.HV = handles.VH;
    elseif ~isempty(handles.HV)&&isempty(handles.VH)
        handles.VH = handles.HV;
    end
    set([handles.CoPolAngEllip handles.CrossPolAngEllip],'Enable','on');
else
    set([handles.CoPolAngEllip handles.CrossPolAngEllip],'Enable','off');
    set(handles.RxAngleEllip,'Value',1);
end

%enable controls
set(handles.HHButton,'Enable','on');
set(handles.HVButton,'Enable','on');
set(handles.VHButton,'Enable','on');
set(handles.VVButton,'Enable','on');
set(handles.CircCoButton,'Enable','on');
set(handles.CircCrossButton,'Enable','on');
set(handles.HHPlusVVButton,'Enable','on');
set(handles.HHMinusVVButton,'Enable','on');
set(handles.Color,'Enable','on');
set(handles.Grayscale,'Enable','on');
set(handles.TxAngleMan,'Enable','on');
set(handles.TxAngleInc,'Enable','on');
set(handles.TxAngleUp,'Enable','on');
set(handles.TxAngleDown,'Enable','on');
set(handles.TxEllipMan,'Enable','on');
set(handles.TxEllipInc,'Enable','on');
set(handles.TxEllipUp,'Enable','on');
set(handles.TxEllipDown,'Enable','on');
set(handles.AngleLockMan,'Enable','on');
set(handles.EllipLockMan,'Enable','on');

set(handles.RunDecomp,'Enable','on');
set(handles.Generate,'Enable','on');
set(handles.SaveDecomp,'Enable','on');
set(handles.SelectPixel,'Enable','on');
set(handles.SelectAOI,'Enable','on');
set(handles.SetRed,'Enable','on');
set(handles.SetGreen,'Enable','on');
set(handles.SetBlue,'Enable','on');

set([handles.TxManPanel handles.RxManPanel],'Visible','on');

%disable controls if only one Tx Pol
if (get(handles.TxHCheck,'Value') == 0)
    %disable appropriate buttons
    set(handles.HHButton,'enable','off');
    set(handles.HVButton,'enable','off');

    set(handles.CircCoButton,'Enable','off');
    set(handles.CircCrossButton,'Enable','off');
    set(handles.HHPlusVVButton,'Enable','off');
    set(handles.HHMinusVVButton,'Enable','off');

    set(handles.AngleLock,'Value',0);
    set(handles.AngleLock,'Enable','off');
    set(handles.EllipLock,'Value',0);
    set(handles.EllipLock,'Enable','off');
    set(handles.TxAngleMan,'Enable','off');
    set(handles.TxAngleInc,'Enable','off');
    set(handles.TxAngleUp,'Enable','off');
    set(handles.TxAngleDown,'Enable','off');
    set(handles.TxEllipMan,'Enable','off');
    set(handles.TxEllipInc,'Enable','off');
    set(handles.TxEllipUp,'Enable','off');
    set(handles.TxEllipDown,'Enable','off');
    set(handles.AngleLockMan,'Enable','off');
    set(handles.EllipLockMan,'Enable','off');
    set(handles.AngleLockMan,'Value',0);
    set(handles.EllipLockMan,'Value',0);  
    set(handles.TxManPanel,'Visible','off');
end
if (get(handles.TxVCheck,'Value') == 0)
    %disable appropriate buttons
    set(handles.VHButton,'enable','off');
    set(handles.VVButton,'enable','off');

    set(handles.CircCoButton,'Enable','off');
    set(handles.CircCrossButton,'Enable','off');
    set(handles.HHPlusVVButton,'Enable','off');
    set(handles.HHMinusVVButton,'Enable','off');

    set(handles.AngleLock,'Value',0);
    set(handles.AngleLock,'Enable','off');
    set(handles.EllipLock,'Value',0);
    set(handles.EllipLock,'Enable','off');
    set(handles.TxAngleMan,'Enable','off');
    set(handles.TxAngleInc,'Enable','off');
    set(handles.TxAngleUp,'Enable','off');
    set(handles.TxAngleDown,'Enable','off');
    set(handles.TxEllipMan,'Enable','off');
    set(handles.TxEllipInc,'Enable','off');
    set(handles.TxEllipUp,'Enable','off');
    set(handles.TxEllipDown,'Enable','off');
    set(handles.AngleLockMan,'Enable','off');
    set(handles.EllipLockMan,'Enable','off');
    set(handles.AngleLockMan,'Value',0);
    set(handles.EllipLockMan,'Value',0);   
    set(handles.TxManPanel,'Visible','off');
end       

%set data to the first co-pol we have
if (handles.HH)
    set(handles.TxAngle,'String',0);
    set(handles.RxAngle,'String',0);
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);
    handles.TxAngleFloat = 0;
    handles.RxAngleFloat = 0;
    handles.TxEllipticityFloat = 0;
    handles.RxEllipticityFloat = 0;
    PlotPolEllipse(handles.TxPlot,0,0);
    PlotPolEllipse(handles.RxPlot,0,0);
    handles = UpdateImage(handles);
else
    set(handles.TxAngle,'String',90);
    set(handles.RxAngle,'String',90);
    set(handles.TxEllipticity,'String',0);
    set(handles.RxEllipticity,'String',0);
    handles.TxAngleFloat = 90;
    handles.RxAngleFloat = 90;
    handles.TxEllipticityFloat = 0;
    handles.RxEllipticityFloat = 0;
    PlotPolEllipse(handles.TxPlot,90,0);
    PlotPolEllipse(handles.RxPlot,90,0);
    handles = UpdateImage(handles);
end

% Update handles structure
guidata(hObject, handles);


function handles = UpdateImage(handles)
% Compute coefficients
[coefs(1,1,handles.HH) coefs(1,1,handles.HV) ...
    coefs(1,1,handles.VH) coefs(1,1,handles.VV)] = ...
    ComputePolCoeff(handles.TxAngleFloat,handles.RxAngleFloat,...
    handles.TxEllipticityFloat,handles.RxEllipticityFloat);
% Apply coefficients
handles.mitm_hand.DataTransformFcn = @(x) sum(bsxfun(@times,x,coefs),3);

if (handles.aviflag == 1)
    %record frame
    temp.TxAngle = handles.TxAngleFloat;
    temp.RxAngle = handles.RxAngleFloat;
    temp.TxEllip = handles.TxEllipticityFloat;
    temp.RxEllip = handles.RxEllipticityFloat;
    temp.Color = get(handles.Color,'Value');
    handles.aviparams = horzcat(handles.aviparams,temp);         
    guidata(gcf, handles); % Update handles structure
end


% --- Executes on button press in AngleLock.
function AngleLock_Callback(hObject, eventdata, handles)
% hObject    handle to AngleLock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AngleLock


% --- Executes on button press in EllipCheck.
function EllipCheck_Callback(hObject, eventdata, handles)
% hObject    handle to EllipCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EllipCheck


% --- Executes on button press in HHButton.
function HHButton_Callback(hObject, eventdata, handles)
% hObject    handle to HHButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',0);
set(handles.RxAngle,'String',0);
set(handles.TxAngleMan,'String',0);
set(handles.RxAngleMan,'String',0);
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
handles.TxAngleFloat = 0;
handles.RxAngleFloat = 0;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
PlotPolEllipse(handles.TxPlot,0,0);
PlotPolEllipse(handles.RxPlot,0,0);
UpdateImage(handles);


% --- Executes on button press in VVButton.
function VVButton_Callback(hObject, eventdata, handles)
% hObject    handle to VVButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',90);
set(handles.RxAngle,'String',90);
set(handles.TxAngleMan,'String',90);
set(handles.RxAngleMan,'String',90);
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
handles.TxAngleFloat = 90;
handles.RxAngleFloat = 90;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
PlotPolEllipse(handles.TxPlot,90,0);
PlotPolEllipse(handles.RxPlot,90,0);
UpdateImage(handles);


% --- Executes on button press in HVButton.
function HVButton_Callback(hObject, eventdata, handles)
% hObject    handle to HVButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',0);
set(handles.RxAngle,'String',90);
set(handles.TxAngleMan,'String',0);
set(handles.RxAngleMan,'String',90);
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
handles.TxAngleFloat = 0;
handles.RxAngleFloat = 90;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
PlotPolEllipse(handles.TxPlot,0,0);
PlotPolEllipse(handles.RxPlot,90,0);
UpdateImage(handles);


% --- Executes on button press in VHButton.
function VHButton_Callback(hObject, eventdata, handles)
% hObject    handle to VHButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',90);
set(handles.RxAngle,'String',0);
set(handles.TxAngleMan,'String',90);
set(handles.RxAngleMan,'String',0);
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
handles.TxAngleFloat = 90;
handles.RxAngleFloat = 0;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
PlotPolEllipse(handles.TxPlot,90,0);
PlotPolEllipse(handles.RxPlot,0,0);
UpdateImage(handles);


% --- Executes on button press in CircCoButton.
function CircCoButton_Callback(hObject, eventdata, handles)
% hObject    handle to CircCoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',0);
set(handles.RxAngle,'String',0);
set(handles.TxAngleMan,'String',0);
set(handles.RxAngleMan,'String',0);
set(handles.TxEllipticity,'String',45);
set(handles.RxEllipticity,'String',45);
set(handles.TxEllipMan,'String',45);
set(handles.RxEllipMan,'String',45);
handles.TxAngleFloat = 0;
handles.RxAngleFloat = 0;
handles.TxEllipticityFloat = 45;
handles.RxEllipticityFloat = 45;
PlotPolEllipse(handles.TxPlot,0,45);
PlotPolEllipse(handles.RxPlot,0,45);
UpdateImage(handles);


% --- Executes on button press in CircCrossButton.
function CircCrossButton_Callback(hObject, eventdata, handles)
% hObject    handle to CircCrossButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',0);
set(handles.RxAngle,'String',0);
set(handles.TxAngleMan,'String',0);
set(handles.RxAngleMan,'String',0);
set(handles.TxEllipticity,'String',45);
set(handles.RxEllipticity,'String',-45);
set(handles.TxEllipMan,'String',45);
set(handles.RxEllipMan,'String',-45);
handles.TxAngleFloat = 0;
handles.RxAngleFloat = 0;
handles.TxEllipticityFloat = 45;
handles.RxEllipticityFloat = -45;
PlotPolEllipse(handles.TxPlot,0,45);
PlotPolEllipse(handles.RxPlot,0,-45);
UpdateImage(handles);


% --- Executes on button press in HHPlusVVButton.
function HHPlusVVButton_Callback(hObject, eventdata, handles)
% hObject    handle to HHPlusVVButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',135);
set(handles.RxAngle,'String',135);
set(handles.TxAngleMan,'String',135);
set(handles.RxAngleMan,'String',135);
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
handles.TxAngleFloat = 135;
handles.RxAngleFloat = 135;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
PlotPolEllipse(handles.TxPlot,135,0);
PlotPolEllipse(handles.RxPlot,135,0);
UpdateImage(handles);


% --- Executes on button press in HHMinusVVButton.
function HHMinusVVButton_Callback(hObject, eventdata, handles)
% hObject    handle to HHMinusVVButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TxAngle,'String',45);
set(handles.RxAngle,'String',135);
set(handles.TxAngleMan,'String',45);
set(handles.RxAngleMan,'String',135);
set(handles.TxEllipticity,'String',0);
set(handles.RxEllipticity,'String',0);
set(handles.TxEllipMan,'String',0);
set(handles.RxEllipMan,'String',0);
handles.TxAngleFloat = 45;
handles.RxAngleFloat = 135;
handles.TxEllipticityFloat = 0;
handles.RxEllipticityFloat = 0;
PlotPolEllipse(handles.TxPlot,45,0);
PlotPolEllipse(handles.RxPlot,135,0);
UpdateImage(handles);


% --- Executes on button press in AngleControl.
function AngleControl_Callback(hObject, eventdata, handles)
% hObject    handle to AngleControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AngleControl


% --- Executes on button press in EllipControl.
function EllipControl_Callback(hObject, eventdata, handles)
% hObject    handle to EllipControl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EllipControl


function Colormap_Callback(hObject, eventdata)

handles = guidata(hObject);
if get(handles.Color,'Value')
    colormap(handles.mitm_hand.AxesHandle,jet);
else
    colormap(handles.mitm_hand.AxesHandle,gray);
end


% --- Executes on button press in RxHCheck.
function RxHCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RxHCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RxVCheck.
function RxVCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RxVCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RxVCheck


% --- Executes on button press in TxHCheck.
function TxHCheck_Callback(hObject, eventdata, handles)
% hObject    handle to TxHCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.TxHCheck,'Value');
VVal = get(handles.TxVCheck,'Value');
if (val == 0)    
    if (VVal == 0)
        msgbox('Must Have at Least One Tx Polarization!!!');
        set(handles.TxHCheck,'Value',1);
        return;
    else
        %set TxPol to V
        set(handles.TxAngle,'String',90);
        set(handles.TxEllipticity,'String',0);        
        handles.TxAngleFloat = 90;        
        handles.TxEllipticityFloat = 0;        
        PlotPolEllipse(handles.TxPlot,90,0);
        UpdateImage(handles);
        
        %disable appropriate buttons
        set(handles.HHButton,'enable','off');
        set(handles.HVButton,'enable','off');
        set(handles.CircCoButton,'Enable','off');
        set(handles.CircCrossButton,'Enable','off');
        set(handles.HHPlusVVButton,'Enable','off');
        set(handles.HHMinusVVButton,'Enable','off');  
        set(handles.TxManPanel,'Visible','off');
        
        set(handles.AngleLock,'Value',0);
        set(handles.AngleLock,'Enable','off');
        set(handles.EllipLock,'Value',0);
        set(handles.EllipLock,'Enable','off');
    end 
else
    set(handles.HHButton,'enable','on');
    set(handles.HVButton,'enable','on');
    if (VVal == 1)
        set(handles.CircCoButton,'Enable','on');
        set(handles.CircCrossButton,'Enable','on');
        set(handles.HHPlusVVButton,'Enable','on');
        set(handles.HHMinusVVButton,'Enable','on');   
        set(handles.AngleLock,'Enable','on');        
        set(handles.EllipLock,'Enable','on');
        set(handles.TxManPanel,'Visible','on');
    end        
end


% --- Executes on button press in TxVCheck.
function TxVCheck_Callback(hObject, eventdata, handles)
% hObject    handle to TxVCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.TxVCheck,'Value');
HVal = get(handles.TxHCheck,'Value');
if (val == 0)    
    if (HVal == 0)
        msgbox('Must Have at Least One Tx Polarization!!!');
        set(handles.TxVCheck,'Value',1);
        return;
    else
        %set TxPol to V
        set(handles.TxAngle,'String',0);
        set(handles.TxEllipticity,'String',0);        
        handles.TxAngleFloat = 0;        
        handles.TxEllipticityFloat = 0;        
        PlotPolEllipse(handles.TxPlot,0,0);
        UpdateImage(handles);
        
        %disable appropriate buttons
        set(handles.VHButton,'enable','off');
        set(handles.VVButton,'enable','off');
        set(handles.CircCoButton,'Enable','off');
        set(handles.CircCrossButton,'Enable','off');
        set(handles.HHPlusVVButton,'Enable','off');
        set(handles.HHMinusVVButton,'Enable','off');
        set(handles.TxManPanel,'Visible','off');
        
        set(handles.AngleLock,'Value',0);
        set(handles.AngleLock,'Enable','off');
        set(handles.EllipLock,'Value',0);
        set(handles.EllipLock,'Enable','off');
    end
else
    set(handles.VHButton,'enable','on');
    set(handles.VVButton,'enable','on');
    if (HVal == 1)
        set(handles.CircCoButton,'Enable','on');
        set(handles.CircCrossButton,'Enable','on');
        set(handles.HHPlusVVButton,'Enable','on');
        set(handles.HHMinusVVButton,'Enable','on');
        set(handles.AngleLock,'Enable','on');        
        set(handles.EllipLock,'Enable','on');
        set(handles.TxManPanel,'Visible','on');
    end    
end


% --- Executes on button press in EllipLock.
function EllipLock_Callback(hObject, eventdata, handles)
% hObject    handle to EllipLock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EllipLock


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


% --- Executes on button press in RxAngleUp.
function RxAngleUp_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngleUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RxAng = str2double(get(handles.RxAngleMan,'String')); 
RxEllip = str2double(get(handles.RxEllipMan,'String')); 
RxAngleInc = str2double(get(handles.RxAngleInc,'String')); 
AngleLock = get(handles.AngleLockMan,'Value');

RxAng = RxAng + RxAngleInc;
if (RxAng > 360)
    RxAng = RxAng - 360;
end
handles.RxAngleFloat = RxAng;
set(handles.RxAngleMan,'String',RxAng);
set(handles.RxAngle,'String',round(RxAng));
PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);

if (AngleLock)
    TxAng = str2double(get(handles.TxAngleMan,'String')); 
    TxEllip = str2double(get(handles.TxEllipMan,'String')); 
    TxAng = TxAng + RxAngleInc;
    if (TxAng > 360)
        TxAng = TxAng - 360;
    end
    handles.TxAngleFloat = TxAng;
    set(handles.TxAngleMan,'String',TxAng);
    set(handles.TxAngle,'String',round(TxAng));
    PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);
end

UpdateImage(handles);


% --- Executes on button press in RxAngleDown.
function RxAngleDown_Callback(hObject, eventdata, handles)
% hObject    handle to RxAngleDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxAng = str2double(get(handles.RxAngleMan,'String')); 
RxEllip = str2double(get(handles.RxEllipMan,'String')); 
RxAngleInc = str2double(get(handles.RxAngleInc,'String')); 
AngleLock = get(handles.AngleLockMan,'Value');

RxAng = RxAng - RxAngleInc;
if (RxAng < 0)
    RxAng = RxAng + 360;
end
handles.RxAngleFloat = RxAng;
set(handles.RxAngleMan,'String',RxAng);
set(handles.RxAngle,'String',round(RxAng));
PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);

if (AngleLock)
    TxAng = str2double(get(handles.TxAngleMan,'String')); 
    TxEllip = str2double(get(handles.TxEllipMan,'String')); 
    TxAng = TxAng - RxAngleInc;
    if (TxAng < 0)
        TxAng = TxAng + 360;
    end
    handles.TxAngleFloat = TxAng;
    set(handles.TxAngleMan,'String',TxAng);
    set(handles.TxAngle,'String',round(TxAng));
    PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);
end

UpdateImage(handles);


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


% --- Executes on button press in RxEllipUp.
function RxEllipUp_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxAng = str2double(get(handles.RxAngleMan,'String')); 
RxEllip = str2double(get(handles.RxEllipMan,'String')); 
RxEllipInc = str2double(get(handles.RxEllipInc,'String')); 
EllipLock = get(handles.EllipLockMan,'Value');

RxEllip = RxEllip + RxEllipInc;
if (RxEllip > 45)
    RxEllip = RxEllip - 90;
end
handles.RxEllipticityFloat = RxEllip;
set(handles.RxEllipMan,'String',RxEllip);
set(handles.RxEllipticity,'String',round(RxEllip));
PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);

if (EllipLock)
    TxAng = str2double(get(handles.TxAngleMan,'String')); 
    TxEllip = str2double(get(handles.TxEllipMan,'String')); 
    TxEllip = TxEllip + RxEllipInc;
    if (TxEllip > 45)
        TxEllip = TxEllip - 90;
    end
    handles.TxEllipFloat = TxEllip;
    set(handles.TxEllipMan,'String',TxEllip);
    set(handles.TxEllipticity,'String',round(TxEllip));
    PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);
end

UpdateImage(handles);


% --- Executes on button press in RxEllipDown.
function RxEllipDown_Callback(hObject, eventdata, handles)
% hObject    handle to RxEllipDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RxAng = str2double(get(handles.RxAngleMan,'String')); 
RxEllip = str2double(get(handles.RxEllipMan,'String')); 
RxEllipInc = str2double(get(handles.RxEllipInc,'String')); 
EllipLock = get(handles.EllipLockMan,'Value');

RxEllip = RxEllip - RxEllipInc;
if (RxEllip < -45)
    RxEllip = RxEllip + 90;
end
handles.RxEllipticityFloat = RxEllip;
set(handles.RxEllipMan,'String',RxEllip);
set(handles.RxEllipticity,'String',round(RxEllip));
PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);

if (EllipLock)
    TxAng = str2double(get(handles.TxAngleMan,'String')); 
    TxEllip = str2double(get(handles.TxEllipMan,'String')); 
    TxEllip = TxEllip - RxEllipInc;
    if (TxEllip < -45)
        TxEllip = TxEllip + 90;
    end
    handles.TxEllipFloat = TxEllip;
    set(handles.TxEllipMan,'String',TxEllip);
    set(handles.TxEllipticity,'String',round(TxEllip));
    PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);
end

UpdateImage(handles);


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


% --- Executes on button press in TxAngleUp.
function TxAngleUp_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngleUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxAng = str2double(get(handles.TxAngleMan,'String')); 
TxEllip = str2double(get(handles.TxEllipMan,'String')); 
TxAngleInc = str2double(get(handles.TxAngleInc,'String')); 
AngleLock = get(handles.AngleLockMan,'Value');

TxAng = TxAng + TxAngleInc;
if (TxAng > 360)
    TxAng = TxAng - 360;
end
handles.TxAngleFloat = TxAng;
set(handles.TxAngleMan,'String',TxAng);
set(handles.TxAngle,'String',round(TxAng));
PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);

if (AngleLock)
    RxAng = str2double(get(handles.RxAngleMan,'String')); 
    RxEllip = str2double(get(handles.RxEllipMan,'String')); 
    RxAng = RxAng + TxAngleInc;
    if (RxAng > 360)
        RxAng = RxAng - 360;
    end
    handles.RxAngleFloat = RxAng;
    set(handles.RxAngleMan,'String',RxAng);
    set(handles.RxAngle,'String',round(RxAng));
    PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);
end

UpdateImage(handles);


% --- Executes on button press in TxAngleDown.
function TxAngleDown_Callback(hObject, eventdata, handles)
% hObject    handle to TxAngleDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxAng = str2double(get(handles.TxAngleMan,'String')); 
TxEllip = str2double(get(handles.TxEllipMan,'String')); 
TxAngleInc = str2double(get(handles.TxAngleInc,'String')); 
AngleLock = get(handles.AngleLockMan,'Value');

TxAng = TxAng - TxAngleInc;
if (TxAng < 0)
    TxAng = TxAng + 360;
end
handles.TxAngleFloat = TxAng;
set(handles.TxAngleMan,'String',TxAng);
set(handles.TxAngle,'String',round(TxAng));
PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);

if (AngleLock)
    RxAng = str2double(get(handles.RxAngleMan,'String')); 
    RxEllip = str2double(get(handles.RxEllipMan,'String')); 
    RxAng = RxAng - TxAngleInc;
    if (RxAng < 0)
        RxAng = RxAng + 360;
    end
    handles.RxAngleFloat = RxAng;
    set(handles.RxAngleMan,'String',RxAng);
    set(handles.RxAngle,'String',round(RxAng));
    PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);
end

UpdateImage(handles);


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


% --- Executes on button press in TxEllipUp.
function TxEllipUp_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxAng = str2double(get(handles.TxAngleMan,'String')); 
TxEllip = str2double(get(handles.TxEllipMan,'String')); 
TxEllipInc = str2double(get(handles.TxEllipInc,'String')); 
EllipLock = get(handles.EllipLockMan,'Value');

TxEllip = TxEllip + TxEllipInc;
if (TxEllip > 45)
    TxEllip = TxEllip - 90;
end
handles.TxEllipticityFloat = TxEllip;
set(handles.TxEllipMan,'String',TxEllip);
set(handles.TxEllipticity,'String',round(TxEllip));
PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);

if (EllipLock)
    RxAng = str2double(get(handles.RxAngleMan,'String')); 
    RxEllip = str2double(get(handles.RxEllipMan,'String')); 
    RxEllip = RxEllip + TxEllipInc;
    if (RxEllip > 45)
        RxEllip = RxEllip - 90;
    end
    handles.RxEllipFloat = RxEllip;
    set(handles.RxEllipMan,'String',RxEllip);
    set(handles.RxEllipticity,'String',round(RxEllip));
    PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);
end

UpdateImage(handles);


% --- Executes on button press in TxEllipDown.
function TxEllipDown_Callback(hObject, eventdata, handles)
% hObject    handle to TxEllipDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TxAng = str2double(get(handles.TxAngleMan,'String')); 
TxEllip = str2double(get(handles.TxEllipMan,'String')); 
TxEllipInc = str2double(get(handles.TxEllipInc,'String')); 
EllipLock = get(handles.EllipLockMan,'Value');

TxEllip = TxEllip - TxEllipInc;
if (TxEllip < -45)
    TxEllip = TxEllip + 90;
end
handles.TxEllipticityFloat = TxEllip;
set(handles.TxEllipMan,'String',TxEllip);
set(handles.TxEllipticity,'String',round(TxEllip));
PlotPolEllipse(handles.TxPlot,TxAng,TxEllip);

if (EllipLock)
    RxAng = str2double(get(handles.RxAngleMan,'String')); 
    RxEllip = str2double(get(handles.RxEllipMan,'String')); 
    RxEllip = RxEllip - TxEllipInc;
    if (RxEllip < -45)
        RxEllip = RxEllip + 90;
    end
    handles.RxEllipFloat = RxEllip;
    set(handles.RxEllipMan,'String',RxEllip);
    set(handles.RxEllipticity,'String',round(RxEllip));
    PlotPolEllipse(handles.RxPlot,RxAng,RxEllip);
end

UpdateImage(handles);


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

[fname, path] = uiputfile( {'*.jpg','jpg Files (*.jpg)'},'Save JPG File',pathstr);

setpref('matlab_sar_toolbox','last_used_directory',path); %store path

set(handles.ImageName,'String',strcat(path,fname));
set(handles.SaveImage,'Enable','on');


% --- Executes on button press in SaveImage.
function SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

F = getframe( handles.mitm_hand.AxesHandle );
imwrite(F.cdata,get(handles.ImageName,'String'),'jpg');   
msgbox('Image Saved');


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

[fname, path] = uiputfile( {'*.avi','AVI Files (*.avi)'},'Save AVI File',pathstr);

setpref('matlab_sar_toolbox','last_used_directory',path); %store path

set(handles.MovieName,'String',strcat(path,fname));
set(handles.StartMovie,'Enable','on');


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


% --- Executes on button press in StartMovie.
function StartMovie_Callback(hObject, eventdata, handles)
% hObject    handle to StartMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%we'll just store the Pol values and zoom/contrast/color params and then
%create the movie 
handles.aviparams = [];

%set avi recording flag
handles.aviflag = 1;

set(handles.StartMovie,'Enable','off');
set(handles.StopMovie,'Enable','on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in StopMovie.
function StopMovie_Callback(hObject, eventdata, handles)
% hObject    handle to StopMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.aviflag = 0;

set(handles.StartMovie,'Enable','off');
set(handles.StopMovie,'Enable','off');

%open avi
framerate = str2double(get(handles.FrameRate,'String'));
avifilename = get(handles.MovieName,'String');
mov = avifile( avifilename, 'compression', 'none', 'fps', framerate );

numframes = length(handles.aviparams);
set(handles.TotFrames,'String',numframes);

for i=1:numframes
    set(handles.TxAngle,'String',round(handles.aviparams(i).TxAngle));
    set(handles.RxAngle,'String',round(handles.aviparams(i).RxAngle));
    set(handles.TxAngleMan,'String',handles.aviparams(i).TxAngle);
    set(handles.RxAngleMan,'String',handles.aviparams(i).RxAngle);
    set(handles.TxEllipticity,'String',round(handles.aviparams(i).TxEllip));
    set(handles.RxEllipticity,'String',round(handles.aviparams(i).RxEllip));
    set(handles.TxEllipMan,'String',handles.aviparams(i).TxEllip);
    set(handles.RxEllipMan,'String',handles.aviparams(i).RxEllip);
    handles.TxAngleFloat = handles.aviparams(i).TxAngle;
    handles.RxAngleFloat = handles.aviparams(i).RxAngle;
    handles.TxEllipticityFloat = handles.aviparams(i).TxEllip;
    handles.RxEllipticityFloat = handles.aviparams(i).RxEllip;
    set(handles.Color,'Value',handles.aviparams(i).Color);
    UpdateImage(handles);
    PlotPolEllipse(handles.RxPlot,handles.aviparams(i).RxAngle,...
                   handles.aviparams(i).RxEllip);
    PlotPolEllipse(handles.TxPlot,handles.aviparams(i).TxAngle,...
                   handles.aviparams(i).TxEllip);
    F = getframe( handles.mitm_hand.AxesHandle );
    mov = addframe( mov, F );
    %update status on GUI
    set(handles.CurFrame,'String',i);
end

mov = close(mov);

set(handles.TotFrames,'String','');
set(handles.CurFrame,'String','');

set(handles.StartMovie,'Enable','on');
set(handles.StopMovie,'Enable','off');

% Update handles structure
guidata(hObject, handles);


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


% --- Executes during object creation, after setting all properties.
function StopMovie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StopMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function TotFrames_Callback(hObject, eventdata, handles)
% hObject    handle to TotFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TotFrames as text
%        str2double(get(hObject,'String')) returns contents of TotFrames as a double


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


% --- Executes on button press in PlotResults.
function PlotResults_Callback(hObject, eventdata, handles)
% hObject    handle to PlotResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AngRes = str2double(get(handles.AngRes,'String'));
ylabels = -45:AngRes:45;
xlabels = 0:AngRes:180;
if get(handles.CoPolAngEllip,'Value') % Co-Pol Angle/Ellipse Plot
    TxAngle = xlabels;
    RxAngle = xlabels;
    TxEllip = ylabels;
    RxEllip = ylabels;
    title_str = 'Co-Pol Angle/Ellipticity Plot';
elseif (get(handles.CrossPolAngEllip,'Value')) % Cross-Pol Angle/Ellipse Plot
    TxAngle = xlabels;
    RxAngle = xlabels+90;
    TxEllip = ylabels;
    RxEllip = -ylabels;
    title_str = 'Cross-Pol Angle/Ellipticity Plot';
else % Rx - Angle/Ellip
    TxAngle = ones(1,numel(xlabels))*handles.TxAngleFloat;
    TxEllip = ones(1,numel(ylabels))*handles.TxEllipticityFloat;
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
imagesc(xlabels,ylabels,data);
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


function AngRes_Callback(hObject, eventdata, handles)
% hObject    handle to AngRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AngRes as text
%        str2double(get(hObject,'String')) returns contents of AngRes as a double


% --- Executes during object creation, after setting all properties.
function AngRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AngRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AngRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in DecompList.
function DecompList_Callback(hObject, eventdata, handles)
% hObject    handle to DecompList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DecompList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DecompList


% --- Executes during object creation, after setting all properties.
function DecompList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DecompList (see GCBO)
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
if isempty(handles.HH)
    HHData = zeros(data_size);
else
    HHData = single(handles.mitm_hand.ComplexData(:,:,handles.HH));
end
if isempty(handles.HV)
    HVData = zeros(data_size);
else
    HVData = single(handles.mitm_hand.ComplexData(:,:,handles.HV));
end
if isempty(handles.VH)
    VHData = zeros(data_size);
else
    VHData = single(handles.mitm_hand.ComplexData(:,:,handles.VH));
end
if isempty(handles.VV)
    VVData = zeros(data_size);
else
    VVData = single(handles.mitm_hand.ComplexData(:,:,handles.VV));
end

%run decomposition from list
%get name of function to call
decompindex = get(handles.DecompList,'Value');
decompstrings = get(handles.DecompList,'String');
decomp = decompstrings(decompindex);

%determine if this is a matlab file or a decomp file (decomp file has
%extension)
if (isempty(strfind(decomp{1},'.decomp')))
    RGBData = feval(decomp{1}, HHData, HVData, VHData, VVData);
else
    RGBData = CreateDecomp(decomp{1}, HHData, HVData, VHData, VVData);
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


function RedPol_Callback(hObject, eventdata, handles)
% hObject    handle to RedPol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RedPol as text
%        str2double(get(hObject,'String')) returns contents of RedPol as a double


% --- Executes during object creation, after setting all properties.
function RedPol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedPol (see GCBO)
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
set(handles.RedPol,'String',PolString);
guidata(hObject, handles);


function GreenPol_Callback(hObject, eventdata, handles)
% hObject    handle to GreenPol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GreenPol as text
%        str2double(get(hObject,'String')) returns contents of GreenPol as a double


% --- Executes during object creation, after setting all properties.
function GreenPol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GreenPol (see GCBO)
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
set(handles.GreenPol,'String',PolString);
guidata(hObject, handles);


function BluePol_Callback(hObject, eventdata, handles)
% hObject    handle to BluePol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BluePol as text
%        str2double(get(hObject,'String')) returns contents of BluePol as a double


% --- Executes during object creation, after setting all properties.
function BluePol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BluePol (see GCBO)
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
set(handles.BluePol,'String',PolString);
guidata(hObject, handles);


% --- Executes on button press in Generate.
function Generate_Callback(hObject, eventdata, handles)
% hObject    handle to Generate (see GCBO)
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

set(handles.DecompList,'String',filestring);


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

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////