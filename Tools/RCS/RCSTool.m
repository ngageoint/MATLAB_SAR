function varargout = RCSTool(varargin)
% RCSTool GUI for computing and displaying Slow/Fast time RCS data profiles
%
% Warning: This tool and the functions used within it assume a common
% radiometric scale factor (or 2D polynomial) for all polarimetric
% channels.  It also assumes a common noise estimate across all
% polarimetric channels.  This is true for much data we have seen; however,
% it is certainly not guaranteed to be true.  Data that does not behave
% this way will result in incorrect results, as the tool is currently
% constructed.  To truly fix this, the SICD-like structure returned from
% hg_mitm_viewer would have to be changed, since it only returns a single
% field for these things across all channels.
%
% Authors: Tim Cox, NRL; Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RCSTool_OpeningFcn, ...
                   'gui_OutputFcn',  @RCSTool_OutputFcn, ...
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


% --- Executes just before RCSTool is made visible.
function RCSTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RCSTool (see VARARGIN)

% Choose default command line output for RCSTool
handles.output = hObject;

% Array of handles for each ROI defined
handles.rois = {};
handles.rois_native = {}; % Temporarily store shape info in full image coordinates for pan/zoom
handles.cad = {};          % Added as part of Dr. Zweber CAD overlay - BDK 4/20/2016
handles.cad_native = {};       % Temporarily store shape info in full image coordinates for pan/zoom  %% Added as part of Dr. Zweber CAD overlay - BDK 4/20/2016

% Setup ShapeTable
set(handles.ShapeTable,'ColumnName',{'Name','Shape','Color','Use','RCS','Pred. Noise'});
set(handles.ShapeTable,'ColumnFormat',{'char','char','char',[],'char','char'});  
set(handles.ShapeTable,'ColumnEditable',[true false false true false false]);  
set(handles.ShapeTable,'ColumnWidth',{75 75 62 40 300 75});
set(handles.ShapeTable,'Data',cell(0,6));

% Setup image display
handles.mitm_hand = hg_mitm_viewer(handles.image);
handles.mitm_hand.PreChangeViewFcn = @() saveShapesNativeCoords(hObject);
handles.mitm_hand.PostChangeViewFcn = @() restoreShapesLocalCoords(hObject);

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[]);
p.addParamValue('segment',1);
p.parse(varargin{:});

guidata(hObject, handles); % Update handles structure before LoadImage

if ~isempty(p.Results.filename)
    LoadImage(p.Results.filename,handles,p.Results.aoi,p.Results.segment);     
end


% --- Outputs from this function are returned to the command line.
function varargout = RCSTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


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

%get complex image
[filename, pathstr] = uigetfile(sar_file_extensions('complex'),'Select Input File',pathstr,'MultiSelect', 'on');   
fullfilename = {};
if iscell(filename)
    for j=1:length(filename)
        fullfilename{j}=[pathstr filename{j}];
    end
elseif filename
    fullfilename{1}=[pathstr filename]; % Assure consistent input format
else
    return; % Cancel was selected
end
handles = DeleteAll_Callback(hObject, eventdata, handles);
LoadImage(fullfilename,handles,[],1);

setpref('matlab_sar_toolbox','last_used_directory',pathstr); %store path

% Update handles structure
guidata(hObject, handles);

function LoadImage(filenames,handles,AOI,segment)

if iscell(filenames)
    set(handles.filename,'String',filenames{1});
else
    set(handles.filename,'String',filenames);
end
handles.mitm_hand.close();
handles.mitm_hand.Frame = segment;
if isempty(AOI) % Default is to show entire image
    handles.mitm_hand.openFile(filenames);
else % AOI was passed in
    handles.mitm_hand.openFile(filenames, true);
    pixels_available = floor(getpixelposition(handles.mitm_hand.AxesHandle))-1;
    handles.mitm_hand.setView('CenterPos',AOI(1:2)+(AOI(3:4)/2),...
        'Zoom',max(ceil(AOI(3:4)./pixels_available(3:4))));
end
meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
if isfield(meta,'Radiometric') && isfield(meta.Radiometric,'RCSSFPoly')
    set(handles.useCalConstVal,'String',num2str(10*log10(meta.Radiometric.RCSSFPoly(1))));
else
    set(handles.useCalConstVal,'String','0');
end

% Need to disable some controls if we can't convert pixels to lat/lon
if isempty(point_slant_to_ground([1; 1], ...
        handles.mitm_hand.Metadata{handles.mitm_hand.Frame}))
    set([handles.SaveKML handles.LoadShapes handles.SaveShapes],'enable','off');
else
    set([handles.SaveKML handles.LoadShapes handles.SaveShapes],'enable','on');
end

MetaIcon(handles.mitm_hand.Metadata{handles.mitm_hand.Frame},...
    'handle',handles.metaicon);


% --- Executes on button press in ComputeRCS.
function ComputeRCS_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeRCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~isempty(get(handles.ShapeTable,'Data')))   % Error check to see if data is available in the shape table - BDK 3/29/2016
  if(~isfield(handles,'fig_H'))  % Add field to handles if it doesn't exist - BDK 3/22/2016
    handles.fig_H = [];
  end

  [SlowData,FastData,SlowXData,FastXData,RangeData,RcsFlags] = ComputeRCSData(handles);
  h = PlotRCSData(handles,SlowData,FastData,SlowXData,FastXData,[],RcsFlags);  % Added 'h =' to  line - BDK 3/22/2016
  handles.fig_H = h;  % Assign active figure handles to structure. - BDK 3/22/2016

  % Update handles structure
  guidata(hObject, handles);   % Added to keep track of RCS plot figure handles - BDK 3/22/2016
end

% Compute all data for display/saving to file
function [SlowData,FastData,SlowXData,FastXData,RangeData,RcsFlags] = ComputeRCSData(handles)  % Added Dr. Zweber modification to output variable list - BDK 4/20/2016
% This function computes total RCS, slow/fast-time profiles, and range
% profiles.  Often we only need one of the three, but we call this function
% and compute all three every time we need only one.  This is not an
% efficient way to do things, but we are too lazy to clean up at the
% moment, since all 3 are currently so tightly entangled with each other.
% In hindsight, these computations should have been written more modularly.

meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};

% Check which kind of auto-focus was applied
if isfield(meta,'ImageFormation')&&isfield(meta.ImageFormation,'AzAutofocus')&&...
        strcmpi(meta.ImageFormation.AzAutofocus,'sv')
    warning('RCSTOOL:SV_AUTOFOCUS','Spatially variant autofocus may result in unexpected results!');
end

% For absolute RCS
[cal_sf, RcsFlags.calfactorset] = get_cal_factor(handles);

% Compute fast time axis labels.  Unlike slow-time axis, fast time axis are
% same throughout image, so we can compute once outside of loop.
[fast_range, RcsFlags.axis_labels.x.fast] = fast_axis(meta);

% Compute labels for measurement units
measure_type_strings = get(handles.measuretype,'String');
measure_type = measure_type_strings{get(handles.measuretype,'Value')};
RcsFlags.axis_labels.y = measure_label(measure_type, ...
    get(handles.dBCheck, 'Value'), RcsFlags.calfactorset);

% Get polarizations.  For a single sicd, the spec says this should always
% be a 3-character string, but viewer sometimes makes this a cell array to
% handle combined datasets.
pols = meta.ImageFormation.TxRcvPolarizationProc;
if ~iscell(pols), pols = {pols}; end

% Process shape table data
shape_table_data = get(handles.ShapeTable,'data');
use_indices = find([shape_table_data{:,4}]); % Use only the shapes checked "Use"
[SlowData, FastData, SlowXData, FastXData, RangeData, totalRCS, nz_data_start, nz_data_stop]=...
    deal(cell(1,numel(use_indices))); % Initialize output data  %%% Added Dr. Zweber modifications of RangeData to the call - BDK 4/20/2016
if isempty(use_indices), return; end
for i=1:numel(use_indices)
    % Getting shape positions
    switch shape_table_data{use_indices(i),2} % Shape type
        case 'Rectangle'
            pos = handles.rois{use_indices(i)}.getPosition();
            pos = [pos(1), pos(2); pos(1)+pos(3), pos(2);...
                pos(1)+pos(3), pos(2)+pos(4); pos(1), pos(2)+pos(4)];
        case 'Polygon'
            pos = handles.rois{use_indices(i)}.getPosition();
        case 'Ellipse'
            pos = handles.rois{use_indices(i)}.getVertices();
    end
    for j = 1:size(pos,1)
        pos(j,:) = handles.mitm_hand.axescoords2native(pos(j,:));
    end
    % Get complex data that covers all vertices
    left = max(1,floor(min(pos(:,1))));
    top = max(1,floor(min(pos(:,2))));
    right = min(double(meta.ImageData.NumCols), ceil(max(pos(:,1))));
    bottom = min(double(meta.ImageData.NumRows), ceil(max(pos(:,2))));
    complex_data = handles.mitm_hand.readerobj{handles.mitm_hand.Frame}.read_chip...
        ([left right],[top bottom]);
    complex_data = single(complex_data); % Required for doing math on data
    mask = poly2mask(pos(:,2)-top+1, pos(:,1)-left+1, ...
        right-left+1, bottom-top+1);

    % Compute radiometric power
    OSR = [1./(meta.Grid.Col.SS * meta.Grid.Col.ImpRespBW), ... % same as AzPad 
           1./(meta.Grid.Row.SS * meta.Grid.Row.ImpRespBW)];
    if ~isscalar(cal_sf) % Polynomial representing a spatially varying RCS factor
        % Convert to a per-pixel scale factor if cal_sf is polynomial
        cal_sf_per_pixel = sicd_polyval2d(cal_sf, left:right, top:bottom, meta);
    else
        cal_sf_per_pixel = cal_sf;
    end
    totalRCS{i} = RCS_Compute(...  % Do actual RCS calculations here
        complex_data, OSR, cal_sf_per_pixel, mask);
    [SlowData{i}, FastData{i}] = RCS_ST_FT(...  % Slow/fast-time curves
        complex_data, meta.SCPCOA.SideOfTrack, OSR, cal_sf_per_pixel, mask, meta);
    % Some data requires preconditioning.  For these cases, recompute.
    % Deweighting that might be applied here does not affect absolute RCS,
    % since the average RCS across aperture in recompute_with_fixed_data
    % will be adjusted to match total RCS computed with original data.
    if ~is_normalized_sicd(meta,1)
        SlowData{i} = recompute_with_fixed_data(1);
    end
    if ~is_normalized_sicd(meta,2) % Data needs some preconditioning
        FastData{i} = recompute_with_fixed_data(2);
    end
    RangeData{i} = RCS_Range_Profile(... % Added Dr. Zweber RangeData to call - BDK 4/20/2016
        complex_data, OSR, cal_sf_per_pixel, mask);

    % Determine non-zero data region and guard band sizes
    fftsize=[length(SlowData{i}) length(FastData{i})];
    nz_data_points = fix(fftsize./OSR); 
    nz_data_points = nz_data_points + mod(nz_data_points,2); % force to be even
    guard_band = (fftsize - nz_data_points) / 2;
    nz_data_start{i} = guard_band + 1;
    nz_data_stop{i}  = guard_band + nz_data_points;

    % Scale factors to account for weighting
    if isfield(meta.Grid.Row, 'WgtFunct')
        rng_wght_f=mean(meta.Grid.Row.WgtFunct.^2) / ...
            (mean(meta.Grid.Row.WgtFunct).^2);
    else  % If no weight in metadata, assume uniform weighting
        rng_wght_f=1.0;
    end
    if isfield(meta.Grid.Col, 'WgtFunct')
        az_wght_f=mean(meta.Grid.Col.WgtFunct.^2) / ...
            (mean(meta.Grid.Col.WgtFunct).^2);
    else  % If no weight in metadata, assume uniform weighting
        az_wght_f=1.0;
    end
    area_sp = sum(mask(:)) * meta.Grid.Row.SS * meta.Grid.Col.SS * ...
        (rng_wght_f * az_wght_f); % We wrap the weighting factor into this for convenience
    switch upper(measure_type)
        case 'RCS'
            % No weighting compensation here since totalRCS was computed
            % with weighting, although slow/fast-time curves removed
            % weighting (but were then scaled apprpriately).
            scale_factor = 1;
        case 'SIGMA-0'
            % Normalize by ground area if sigma-0 requested
            scale_factor = cosd(meta.SCPCOA.SlopeAng) / area_sp;
        case 'BETA-0'
            scale_factor = 1 / area_sp;
        case 'GAMMA-0'
            scale_factor = cosd(meta.SCPCOA.SlopeAng) / ...
                (sind(meta.SCPCOA.GrazeAng) * area_sp);
        case 'AVG. PIXEL POWER'  % Raw pixel power
            scale_factor = prod(OSR) / sum(mask(:));
    end
    SlowData{i} = SlowData{i}*scale_factor;
    FastData{i} = FastData{i}*scale_factor;
    totalRCS{i} = totalRCS{i}*scale_factor;
    RangeData{i} = RangeData{i}*scale_factor;   %%% Added Dr. Zweber range data calculation - BDK 4/20/2016
    
    % Axis data
    st_unit_strings = get(handles.slowtimeunits,'String');
    [slow_range, RcsFlags.axis_labels.x.slow] = slow_axis(meta, ...
        st_unit_strings{get(handles.slowtimeunits,'Value')}, ...
        round([(left+right)/2 (top+bottom)/2]), ... % Center of AOI
        str2double(get(handles.target_az,'string'))); % Can vary by location
    SlowXData{i} = linspace(slow_range(1),slow_range(2),length(SlowData{i}));
    FastXData{i} = linspace(fast_range(1),fast_range(2),length(FastData{i}));

    % Compute predicted noise
    if isfield(meta,'Radiometric') && ...
            isfield(meta.Radiometric, 'NoiseLevel') && ...
            isfield(meta.Radiometric.NoiseLevel, 'NoisePoly') && ...
            isfield(meta.Radiometric.NoiseLevel, 'NoiseLevelType') && ...
            strcmpi(meta.Radiometric.NoiseLevel.NoiseLevelType, 'ABSOLUTE') && ...
            any(isfield(meta.Radiometric, ...
            {'RCSSFPoly','SigmaZeroSFPoly','BetaZeroSFPoly','GammaZeroSFPoly'}))
        meta = derived_sicd_fields(meta); % Assure all radiometric fields exist
        switch upper(measure_type)
            case 'RCS'
                sfpoly = meta.Radiometric.RCSSFPoly * sum(mask(:));
            case 'SIGMA-0'
                sfpoly = meta.Radiometric.SigmaZeroSFPoly;
            case 'BETA-0'
                sfpoly = meta.Radiometric.BetaZeroSFPoly;
            case 'GAMMA-0'
                sfpoly = meta.Radiometric.GammaZeroSFPoly;
            case 'AVG. PIXEL POWER'
                sfpoly = 1;
        end
        pred_noise = (10^(sicd_polyval2d(meta.Radiometric.NoiseLevel.NoisePoly, ...
            mean([left right]), mean([top bottom]), meta)/10)) * ... % NoisePoly in dB scale, convert to linear
            sicd_polyval2d(sfpoly, ...
            mean([left right]), mean([top bottom]), meta); % Already in linear scale
        if get(handles.dBCheck,'Value') % Convert to dB
            pred_noise = 10*log10(pred_noise);
        end
        shape_table_data{use_indices(i),6} = num2str(pred_noise);
    end
end
        
% For dB relative RCS data, we scale by total_min.  We wish to maintain
% relative RCS within all ROIs in this loop, so we use a single scaling
% factor for all ROIs.
if RcsFlags.calfactorset
    total_min = ones(size(totalRCS{1})); % Data is calibrated.  Don't do relative scaling.
else
    % A few options:
    % total_min = 1; % No scaling.  User could also do this by forcing cal constant to zero.
    % total_min = min(cellfun(@min,[SlowData FastData])); % Lowest value in all data (likely in zeropad area) is zero dB
    total_min = min(cell2mat(totalRCS),[],2); % Lowest RCS ROI is zero dB
end

% Compute the avg, min and max RCS over the non-zeropad data
for i=1:numel(SlowData)
    if get(handles.dBCheck,'Value') % Convert to dB
        SlowData{i} = 10*log10(bsxfun(@rdivide,SlowData{i},total_min.'));
        FastData{i} = 10*log10(bsxfun(@rdivide,FastData{i},total_min.'));
        RangeData{i} = 10*log10(bsxfun(@rdivide,RangeData{i},total_min.'));  % Added Dr. Zweber modification for range  - BDK 4/20/2016

        RcsFlags.stats(i).slow.avg = 10*log10(totalRCS{i}./total_min);
        RcsFlags.stats(i).fast.avg = 10*log10(totalRCS{i}./total_min);
    else
        RcsFlags.stats(i).slow.avg = totalRCS{i};
        RcsFlags.stats(i).fast.avg = totalRCS{i};
    end
    RcsFlags.stats(i).slow.max = max(SlowData{i}(nz_data_start{i}(1):nz_data_stop{i}(1))); 
    RcsFlags.stats(i).slow.min = min(SlowData{i}(nz_data_start{i}(1):nz_data_stop{i}(1))); 
    RcsFlags.stats(i).fast.max = max(FastData{i}(nz_data_start{i}(2):nz_data_stop{i}(2))); 
    RcsFlags.stats(i).fast.min = min(FastData{i}(nz_data_start{i}(2):nz_data_stop{i}(2))); 

    % Compute guard band edges
    RcsFlags.stats(i).slow.guard_band_low=SlowXData{i}(nz_data_start{i}(1));
    RcsFlags.stats(i).slow.guard_band_hi =SlowXData{i}(nz_data_stop{i}(1));

    RcsFlags.stats(i).fast.guard_band_low=FastXData{i}(nz_data_start{i}(2));
    RcsFlags.stats(i).fast.guard_band_hi =FastXData{i}(nz_data_stop{i}(2));
end

col_names = get(handles.ShapeTable,'ColumnName');
col_names{5} = RcsFlags.axis_labels.y;
set(handles.ShapeTable,'ColumnName',col_names);

for i = 1:numel(use_indices)                        % Loop over each annotation - BDK added on 3/21/2016
    rcs_strs = {};
    for k = 1:numel(RcsFlags.stats(i).slow.avg)
        rcs_strs{end+1} = [num2str(RcsFlags.stats(i).slow.avg(k)) ' ' pols{k}];
    end
    shape_table_data{use_indices(i),5} = strjoin(rcs_strs, ', ');
end
set(handles.ShapeTable,'Data',shape_table_data);    % Put values in shape table. - BDK added on 3/21/2016


% Get calibration scale factor from GUI
function [cal_sf, cal_bool] = get_cal_factor(handles)
    meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
    measure_type_strings = get(handles.measuretype,'String');
    if strcmpi(measure_type_strings{get(handles.measuretype,'Value')},'AVG. PIXEL POWER')
        cal_sf = 1;  % Measure pure pixel power
    elseif get(handles.userCalConst,'Value')
        % User wants to set their own calibration factor, perhaps overriding
        % the one provided with the data.
        calconstdB = str2double(get(handles.useCalConstVal,'string'));
        cal_sf   = 10.^(calconstdB/10.0);
    elseif isfield(meta,'Radiometric') && any(isfield(meta.Radiometric, ...
            {'RCSSFPoly','SigmaZeroSFPoly','BetaZeroSFPoly','GammaZeroSFPoly'}))
        % Use calibration information provided with data
        if ~isfield(meta.Radiometric,'RCSSFPoly')
            meta = derived_sicd_fields(meta); % Derive RCSSFPoly from another field
        end
        cal_sf = meta.Radiometric.RCSSFPoly;
    else
        % Throw warning if no calibration factor is provided
        warning('RCSTOOL:NO_CAL_INFO','No calibration constant provided. RCS results will be relative!');
    end
    cal_bool=exist('cal_sf','var'); % Is calibration factor set?
    if ~cal_bool, cal_sf = 1; end; % Allows same equations to be used


% Compute information we need for annotating slow time axis based on data
% generation options
function [range, label] = slow_axis(meta, units, center_rowcol, relative_az)
try % In case required metadata is not populated
    AzPad = 1/(meta.Grid.Col.ImpRespBW*meta.Grid.Col.SS);
    if strcmpi(units,'Collect Time')
        % Angular aperture spanned by this complex image (without
        % zeropad).  Theta is expressed as a ratio, rather than an angle.
        Theta = meta.Grid.Col.ImpRespBW/meta.Grid.Row.KCtr;
        % This method of estimating effective duration is generic for
        % nearly all types of SAR data, but really only a rough estimate.
        % We compute the effective duration at scene center point only. The
        % effective duration can actually vary slightly across the collect,
        % and we could compute instanteous effective duration at center of
        % aperture for this point at center_rowcol.  However, computing
        % time (as opposed to angle) from complex data is already fairly
        % imprecise, so we don't bother with this additional computation.
        V = norm([meta.SCPCOA.ARPVel.X meta.SCPCOA.ARPVel.Y meta.SCPCOA.ARPVel.Z]);
        EffectiveDuration = Theta * (meta.SCPCOA.SlantRange/V)/...
            sind(meta.SCPCOA.DopplerConeAng);
        % A more accurate way to get total time for PFA spotlight data:
        if isfield(meta,'ImageFormation') && ...
                all(isfield(meta.ImageFormation,{'ImageFormAlgo','TStartProc','TEndProc'})) &&...
                strcmpi(meta.ImageFormation.ImageFormAlgo,'PFA') &&...
                isfield(meta,'CollectionInfo') &&...
                isfield(meta.CollectionInfo,'RadarMode')&&...
                isfield(meta.CollectionInfo.RadarMode,'ModeType')&&...
                strcmpi(meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')               
            EffectiveDuration = meta.ImageFormation.TEndProc - ...
                meta.ImageFormation.TStartProc;
        end
        % Note that using meta.Timeline.CollectDuration is NOT a good way
        % to compute effective duration since an image could be formed from
        % only a portion of a collect and since non-spotlight collects do
        % not illuminative all points in the scene at all times.
        range = [0 EffectiveDuration] + [-1 1]*0.5*EffectiveDuration*(AzPad - 1);
        label = 'Collection Time (sec)';
    elseif strcmpi(units,'Polar Angle')
        % The RCS computation will be ordered with respect to increasing
        % time, not in order of increasing polar angle.  For right-looking
        % collects, polar angle is decreasing with time.
        orientation = [1 -1]; % Normalized bounds of the polar angle
        if meta.SCPCOA.SideOfTrack(1) ~= 'R'
            orientation = - orientation; % Flip for left-looking
        end
        range = atand((orientation./(2*meta.Grid.Col.SS))./meta.Grid.Row.KCtr);
        label = 'Polar Angle (deg)';
    else % Azimuth angle
        st_mid_end = ComputeAz(meta,[0 50 100],center_rowcol); % Start/middle/end degrees
        st_mid_end = unwrap(st_mid_end*pi/180)*180/pi; % Handle cases that cross 0/360
        switch units
            case 'Azimuth Angle'
                relative_az = 0;
                label = 'Azimuth Angle (deg)';
            case 'Aperture Relative' % Relative to center of collect
                relative_az = st_mid_end(2); % Center of aperture
                label = 'Azimuth angle relative to aperture center (deg)';
            case 'Target Relative' % Relative to given target angle
                label = sprintf('Azimuth Angle relative to target at %3.1f (deg)', relative_az);
        end
        % Expand azimuth to include padding
        DeltaAz = diff(st_mid_end([1 end]));
        range = st_mid_end([1 end]) + [-1 1]*0.5*DeltaAz*(AzPad - 1) - relative_az;
    end
catch
    if isfield(meta, 'SCPCOA') && isfield(meta.SCPCOA,'SideOfTrack') && ...
            meta.SCPCOA.SideOfTrack(1) ~= 'R'
        range = [1 0]; % Flip left-looking collects
    else
        range = [0 1];
    end
    label = '';
    warning('RCSTOOL:INSUFFICIENT_METADATA','Insufficient metadata to determine slow-time units.');    
end


% Compute information we need for annotating slow time axis based on data
% generation options
function [range, label] = fast_axis(meta)
    if isfield(meta,'Grid') && isfield(meta.Grid,'Row') && ...
        all(isfield(meta.Grid.Row,{'SS','KCtr'}))
        % Frequency
        BWrg = SPEED_OF_LIGHT/(2*meta.Grid.Row.SS);
        vfrqC = meta.Grid.Row.KCtr*SPEED_OF_LIGHT/2;
        range = (vfrqC + [-1 1]*BWrg/2)/10e8;
        label = 'Frequency (GHz)';
    else % Resort to unitless axis
        range = [0 1];
        label = '';
        warning('RCSTOOL:INSUFFICIENT_METADATA','Insufficient metadata to determine receive frequencies.');
    end

function label = measure_label(measure_type, db, calibrated)
    label = measure_type;
    if strncmpi(measure_type,'RCS',3)
        if calibrated
            if db
                label = [label ' (dBsm)']; 
            else
                label = [label ' m^2']; 
            end
        else
            if db
                label = [label ' (dB)'];
            end
        end
    else % Sigma-0, Beta-0, Gamma-0
        if db
            label = [label ' (dB)']; 
        end
    end
    if ~calibrated
        label = ['Relative ' label];
    end


% Helper function to precondition data and recompute RCS plots
function new_rcs = recompute_with_fixed_data(dim)
    % This is UGLY!.  Mixing namespaces here.  But we don't want to
    % explicitly pass all the variables necessary, but also don't want to
    % have redundant code with this being called twice.  And since GUIDE
    % doesn't produce functions with ENDs, we can't make this a nested
    % function to share the namespace either.  So heavy use of
    % evalin('caller',...) here.
    if dim == 1
        grid_struct = evalin('caller','meta.Grid.Col');
    else
        grid_struct = evalin('caller','meta.Grid.Row');
    end
    
    [ DeltaKCOAPoly, az_coords_m, rg_coords_m ] = deskewparams(...
        evalin('caller','meta'), dim);
    az_coords_m = az_coords_m(evalin('caller','left:right'));
    rg_coords_m = rg_coords_m(evalin('caller','top:bottom'));
    try
        weight_fun = sicdweight2fun(grid_struct);
    catch % Unable to determine weighting function from metadata.
        % Estimate weighting from complex data instead
        weight_fun = estimate_weighting_file(...
            evalin('caller','handles.mitm_hand.readerobj{handles.mitm_hand.Frame}'), dim);
    end
    if isfield(grid_struct,'Sgn')
        fft_sign = grid_struct.Sgn;
    else
        fft_sign = -1;
    end
    complex_data_norm = evalin('caller','complex_data');
    for i = 1:size(complex_data_norm,3) % Handle polarimetric data
        complex_data_norm(:,:,i) = normalize_complex_mem(complex_data_norm(:,:,i),...
            DeltaKCOAPoly, az_coords_m, rg_coords_m, ... % Deskew parameters
            weight_fun, evalin('caller',['OSR(' num2str(dim) ')']), ... % Deweighting parameters
            fft_sign, dim);
    end
    rcs_norm = RCS_Compute(complex_data_norm, evalin('caller','OSR'), ...
        evalin('caller','cal_sf_per_pixel'), evalin('caller','mask')); % Recompute with fixed data
    [rcsdata{1}, rcsdata{2}] = RCS_ST_FT(complex_data_norm, ...
        evalin('caller','meta.SCPCOA.SideOfTrack'), evalin('caller','OSR'), ...
        evalin('caller','cal_sf_per_pixel'), evalin('caller','mask')); % Recompute with fixed data
    new_rcs = rcsdata{dim};
    new_rcs = bsxfun(@times,new_rcs,(evalin('caller','totalRCS{i}')./rcs_norm).');


% Display plots in MATLAB figures
function h = PlotRCSData(handles,SlowData,FastData,SlowXData,FastXData,ShapeNum,RcsFlags)
% Get shape table data
data = get(handles.ShapeTable,'data');
if exist('ShapeNum','var') && ~isempty(ShapeNum)
    indices = ShapeNum;
else
    indices = find([data{:,4}]); % Use only the shapes checked "Use"
end

% Get current RCS plot figure handles - BDK on 3/23/2016 -- Modified to only use good figure handles 6/8/2016
FigHand = handles.fig_H(ishandle(handles.fig_H));

for i=1:numel(indices)
    h = figure;
    FigHand = [FigHand h];   % Build figure handle vector - BDK 3/22/2016
    set(h,'Name',data{indices(i),1});   % Added annotation name to Figure - BDK on 3/22/2016
    % Slow time
    ax_hand = subplot(2,1,1);
    line_hand = plot(ax_hand,SlowXData{indices(i)},SlowData{indices(i)});
    if isvector(SlowData{indices(i)}) % Doesn't work for polarimetric data
        line([SlowXData{indices(i)}(1),SlowXData{indices(i)}(end)],...
            [RcsFlags.stats(indices(i)).slow.avg,RcsFlags.stats(indices(i)).slow.avg],...
            'LineStyle','--','Color','r'); 
        line([SlowXData{indices(i)}(1),SlowXData{indices(i)}(end)],...
            [RcsFlags.stats(indices(i)).slow.max,RcsFlags.stats(indices(i)).slow.max],...
            'LineStyle',':','Color','r');
        line([SlowXData{indices(i)}(1),SlowXData{indices(i)}(end)],...
            [RcsFlags.stats(indices(i)).slow.min,RcsFlags.stats(indices(i)).slow.min],...
            'LineStyle',':','Color','r');
    else % Polarimetric data
        for j = 1:numel(line_hand)
            line([SlowXData{indices(i)}(1),SlowXData{indices(i)}(end)],...
                [RcsFlags.stats(indices(i)).slow.avg(j),RcsFlags.stats(indices(i)).slow.avg(j)],...
                'LineStyle','--','Color',get(line_hand(j),'Color')); 
        end
        legend(ax_hand,handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageFormation.TxRcvPolarizationProc{:})
    end
    grid on; 
    if isempty(RcsFlags.axis_labels.x.slow)
        set(ax_hand,'XTick',[]); % Units unknown
    else
        xlabel(ax_hand,RcsFlags.axis_labels.x.slow);
    end
    ylabel(ax_hand,RcsFlags.axis_labels.y);
    title(ax_hand, ['Slow Time Response: ' data{indices(i),1}]);
    if get(handles.slowtimeunits,'Value')==1 % Azimuth angle
        addlistener(ax_hand,'XTick','PostSet',... % Constrain angles to be within [0 360]
            @(src,evtdata) set(ax_hand, 'XTickLabel',mod(get(ax_hand, 'XTick'), 360)));
    end
    if get(handles.zeropad,'Value')
        xlim(ax_hand,sort([RcsFlags.stats(indices(i)).slow.guard_band_low RcsFlags.stats(indices(i)).slow.guard_band_hi]));
    else
        line([RcsFlags.stats(indices(i)).slow.guard_band_low,RcsFlags.stats(indices(i)).slow.guard_band_low],...
            [min(SlowData{indices(i)}(~isinf(SlowData{indices(i)}))),max(SlowData{indices(i)}(~isinf(SlowData{indices(i)})))],...
            'LineStyle','--','Color','r'); 
        line([RcsFlags.stats(indices(i)).slow.guard_band_hi,RcsFlags.stats(indices(i)).slow.guard_band_hi],...
            [min(SlowData{indices(i)}(~isinf(SlowData{indices(i)}))),max(SlowData{indices(i)}(~isinf(SlowData{indices(i)})))],...
            'LineStyle','--','Color','r'); 
        xlim(ax_hand, [min(SlowXData{indices(i)}) max(SlowXData{indices(i)})]);
    end

    % Fast time
    ax_hand2 = subplot(2,1,2);
    line_hand2 = plot(ax_hand2,FastXData{indices(i)},FastData{indices(i)});
    if isvector(SlowData{indices(i)}) % Doesn't work for polarimetric data
        line([FastXData{indices(i)}(1),FastXData{indices(i)}(end)],...
            [RcsFlags.stats(indices(i)).fast.avg,RcsFlags.stats(indices(i)).fast.avg],...
            'LineStyle','--','Color','r'); 
        line([FastXData{indices(i)}(1),FastXData{indices(i)}(end)],...
            [RcsFlags.stats(indices(i)).fast.max,RcsFlags.stats(indices(i)).fast.max],...
            'LineStyle',':','Color','r');
        line([FastXData{indices(i)}(1),FastXData{indices(i)}(end)],...
            [RcsFlags.stats(indices(i)).fast.min,RcsFlags.stats(indices(i)).fast.min],...
            'LineStyle',':','Color','r');
    else % Polarimetric data
        for j = 1:numel(line_hand)
            line([FastXData{indices(i)}(1),FastXData{indices(i)}(end)],...
                [RcsFlags.stats(indices(i)).fast.avg(j),RcsFlags.stats(indices(i)).fast.avg(j)],...
                'LineStyle','--','Color',get(line_hand2(j),'Color')); 
        end
        legend(ax_hand2,handles.mitm_hand.Metadata{handles.mitm_hand.Frame}.ImageFormation.TxRcvPolarizationProc{:});
    end
    grid on; 
    if isempty(RcsFlags.axis_labels.x.fast)
        set(ax_hand2,'XTick',[]); % Units unknown
    else
        xlabel(ax_hand2,RcsFlags.axis_labels.x.fast);
    end
    ylabel(ax_hand2,RcsFlags.axis_labels.y);
    title(ax_hand2,['Fast Time Response: ' data{indices(i),1}]);
    if get(handles.zeropad,'Value')
        xlim(ax_hand2,sort([RcsFlags.stats(indices(i)).fast.guard_band_low RcsFlags.stats(indices(i)).fast.guard_band_hi]));
    else
        line([RcsFlags.stats(indices(i)).fast.guard_band_low,RcsFlags.stats(indices(i)).fast.guard_band_low],...
            [min(FastData{indices(i)}(~isinf(FastData{indices(i)}))),max(FastData{indices(i)}(~isinf(FastData{indices(i)})))],...
            'LineStyle','--','Color','r'); 
        line([RcsFlags.stats(indices(i)).fast.guard_band_hi,RcsFlags.stats(indices(i)).fast.guard_band_hi],...
            [min(FastData{indices(i)}(~isinf(FastData{indices(i)}))),max(FastData{indices(i)}(~isinf(FastData{indices(i)})))],...
            'LineStyle','--','Color','r'); 
        xlim(ax_hand2,[min(FastXData{indices(i)}) max(FastXData{indices(i)})]);
    end
  % Removed BDK 6/8/2016 end  
end
h = FigHand;  % Pass back the figure handle vector - BDK 3/22/2016

% --- Executes during object creation, after setting all properties.
function ShapeName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShapeName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DrawShape.
function DrawShape_Callback(hObject, eventdata, handles)
% hObject    handle to DrawShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Colors = {'Red','Yellow','Blue','Green','Magenta','White'};    
newColor = Colors{mod(numel(handles.rois),6)+1};

Name = get(handles.ShapeName,'String');
if isempty(Name)
    Name = sprintf('Shape%d',numel(handles.rois));
end

set(hObject,'Enable','off');
old_visible = get(handles.mitm_hand.movie_toolbar,'Visible'); % ROI functions will make everything in axes visible
if get(handles.Rectangle,'Value')
    H = imrect(handles.mitm_hand.AxesHandle);
    roistring = 'Rectangle';
elseif get(handles.Polygon,'Value')
    H = impoly(handles.mitm_hand.AxesHandle);
    roistring = 'Polygon';
elseif get(handles.Ellipse,'Value')
    H = imellipse(handles.mitm_hand.AxesHandle);
    roistring = 'Ellipse';
end
set(handles.mitm_hand.movie_toolbar,'Visible',old_visible); % reset
set(hObject,'Enable','on');
H.setColor(newColor);
handles.rois{end+1} = H;

%%%%% BDK MOD - 3/17/2016 %%%%%
iName = 0;   % Counter to correct label - Added by BDK on 3/17/2016
tempName = get(handles.ShapeTable,'Data');  % Get the ShapeTable fields - Added by BDK 3/17/2016
NameLogical = sum(strcmp(Name, tempName(:,1)));  % Check if new annotation name matches an existing annotation name - Added by BDK 3/17/2016
while(NameLogical>0)  % Error check for name conflict - Added by BDK 3/17/2016
    if(NameLogical)     % Correct name if conflict exists - Added by BDK 3/17/2016
        Name = sprintf('Shape%d',numel(handles.rois)+iName);
        newColor = Colors{mod(numel(handles.rois),6)+1+iName};  % Modify annotation color if name conflict exists
        H.setColor(newColor);
    end
    NameLogical = sum(strcmp(Name, tempName(:,1)));
end
%%%%% BDK MOD - 3/17/2016 %%%%%

%update table data with new shape
set(handles.ShapeTable,'data',[get(handles.ShapeTable,'data'); ...
    {Name, roistring, newColor, true, ' - ', '-'}]);
% Update quantities in table
H.addNewPositionCallback(@(pos) ComputeRCSData(guidata(hObject)));
ComputeRCSData(handles);

% Reset toolbar buttons
try % Really we should check for each components "toggle-button"ness
    tb_components = handles.mitm_hand.main_toolbar.getComponents();
    for i = 2:5 % Hardcoded because we know where toggle buttons are
        tb_components(i).setSelected(0);
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in DeleteAll.
function handles = DeleteAll_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%% BDK MOD - 3/17/2016 %%%%%%
if(~isfield(handles,'fig_H'))           % Error checking: Make figure handles field if it doesn't currently exist. - BDK 3/23/2016
    handles.fig_H = [];
end
if(~isfield(handles,'fig_H_RP'))        % Error checking: Make figure handles field if it doesn't currently exist. - BDK 4/21/2016
    handles.fig_H_RP = [];
end

handles.fig_H = handles.fig_H(ishandle(handles.fig_H));  % Error check the figure handles vector and remove invalid figure handles. - BDK 6/9/2016
handles.fig_H_RP = handles.fig_H_RP(ishandle(handles.fig_H_RP)); % Error check the figure handles vector and remove invalid figure handles. - BDK 6/9/2016

tableData = get(handles.ShapeTable,'Data');     % Obtain the ShapeTable data - BDK 3/17/2016
SelectedData = [tableData{:,4}];                % Break encapsulation - BDK 3/17/2016
iSelect = find(SelectedData == 1);              % Determine user selected data - BDK 3/17/2016
DumpAnnoName = tableData(SelectedData,1);       % Build vector of annotation names being removed - BDK 3/28/2016

for i = 1:length(iSelect);                      % Loop over annotations to remove - BDK 3/17/2016
    delete(handles.rois{iSelect(i)});             % Remove annotations - BDK 3/17/2016
    figHand = handles.fig_H(strcmp(DumpAnnoName(i),get(handles.fig_H,'Name')));  %% BDK Added handles.fig_H(find...
    figHand_RP = handles.fig_H_RP(strcmp(['Range Profile: ' DumpAnnoName{i}],get(handles.fig_H_RP,'Name')));  %% BDK Added handles.fig_H(find...  - Modified by BDK 4/21/2016
    
    if(~isempty(figHand))
        close(figHand);                                                   % Removed selected RCS plot - BDK 4/18/2016
        handles.fig_H = handles.fig_H(ishandle(handles.fig_H));  % Update figure handle vector after deleting the figure - BDK 4/18/2016
    end
    if(~isempty(figHand_RP))
        close(figHand_RP);                                                        % Removed selected RCS plot - BDK 4/21/2016
        handles.fig_H_RP = handles.fig_H_RP(ishandle(handles.fig_H_RP)); % Update figure handle vector after deleting the figure - BDK 4/21/2016
    end
end
handles.rois = handles.rois(SelectedData);             % Rebuild annotation list kept by user - BDK 3/17/2016
set(handles.ShapeTable,'Data',tableData(~SelectedData,:));       % Update the ShapeTable after removing annotatons - BDK 3/17/2016

% Added if statement for removal of all data because sometimes the tool was locking up when many annotations - BDK 3/17/2016
% were added and deleted.  Using this code to delete everything periodically seemed to stabilize the tool. - BDK 3/17/2016
if all(SelectedData)
    % Loop through all the shape handles and delete them
    handles.rois = {};
    % This shouldn't find anything, but sometimes we get out of sync:
    delete(findobj(handles.mitm_hand.AxesHandle,'Type','hggroup'));
    % Delete all the data in the ShapeTable
    set(handles.ShapeTable,'Data',cell(0,6));
    handles.fig_H = [];    % Remove figure handles since all plots are deleted - BDK 3/23/2016
    handles.fig_H_RP = []; % Remove figure handles since all plots are deleted - BDK 4/21/2016
end
%%%%% END OF BDK MOD - 3/17/2016 %%%%%%

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in LoadShapes.
function LoadShapes_Callback(hObject, eventdata, handles)
% hObject    handle to LoadShapes (see GCBO)
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

%get complex image
[filename, pathstr] = uigetfile( {'*.mat','MATLAB Files (*.mat)'},'Select ROI File',pathstr);
fullfilename = fullfile(pathstr,filename);
try 
    valid_file = ~isempty( who('-file',fullfilename,'shape_struct'));
catch
    valid_file = false;
end
if ~valid_file
    msgbox('Not a valid ROI file.','Warning','warn');
    return;
end
shape_struct = getfield(load(fullfilename, 'shape_struct'),'shape_struct');
data = get(handles.ShapeTable,'data');
for i=1:numel(shape_struct)
    % All ROI files are store in polygons in lat/long space, since we don't
    % know which image plane/pixel space they will be loaded into.
    data{end+1,1} = shape_struct(i).Name;
    data{end,2} = 'Polygon';
    data{end,3} = shape_struct(i).Color;
    data{end,4} = shape_struct(i).Toggle;
    %convert lat lon to pixel space
    pos = point_ground_to_slant(shape_struct(i).Vertices,...
        handles.mitm_hand.Metadata{handles.mitm_hand.Frame});
    pos = pos([2 1],:).'; % convert from [row; column] to [column,row]
    for j = 1:size(pos,1)
        pos(j,:) = handles.mitm_hand.native2axescoords(pos(j,:));
    end
    handles.rois{end+1} = impoly(handles.mitm_hand.AxesHandle, pos);
    handles.rois{end}.setColor(shape_struct(i).Color);
    
%%% Begining of CAD overlay code by Peter Zweber - BDK 4/20/2016 %%%
    if ~isempty(regexpi(shape_struct(i).Name,'[l][a][y][o][v][e][r]'))
        % Append handles.cad struct (struct or class?)
        handles.cad{end+1}.index = size(data,1);

        %PKZ (20160624) 
        %The number of fields in shape_struct created by
        %LayoverLeopard, effective 20160624, is 9.
        %Number of fields in shape_struct created by RCSTool SaveShapes is 4.
        %The extra 5 fields from LayoverLeopard deal with the stored CAD model.
        %The '(numel(fieldnames(shape_struct)) == 9)' prevents RCSTool from
        %erroring out when handling a non-LayoverLeopard shape_struct.
        
        if (numel(fieldnames(shape_struct)) == 9)
            %convert lat lon to pixel space
            [layoverVertices] = point_ground_to_slant(shape_struct(i).CADLLAVertices,handles.mitm_hand.Metadata{handles.mitm_hand.Frame});
            layoverVertices = layoverVertices([2 1],:).';   % convert from [row; column] to [column,row]
            layoverVertices = handles.mitm_hand.native2axescoords(layoverVertices); % convert from image coords to axes coords

            % Load CAD model from shape structure - via PlotRangeProfile
            %handles.cad{end}.vertices = shape_struct(i).CADFullImageVertices;
            handles.cad{end}.vertices = layoverVertices;
            handles.cad{end}.faces = shape_struct(i).CADFaces;
            handles.cad{end}.faceAlpha = shape_struct(i).FaceAlpha;
            handles.cad{end}.color = shape_struct(i).Color;
            handles.cad{end}.edgeColor = shape_struct(i).EdgeColor;
            %handles.cad{end}.visibility = 'On';

            %Create the layover patch
            layoverH = patch('vertices',layoverVertices,'faces',shape_struct(i).CADFaces,'facecolor',shape_struct(i).Color,...
                'edgecolor',shape_struct(i).EdgeColor,'FaceAlpha',shape_struct(i).FaceAlpha);
            set(layoverH, 'Visible', 'off'); %PKZ: Can also add to patch declaration
            handles.cad{end}.patch = layoverH;
        end % END of 'numel(fieldnames(shape_struct)) == 9' IF
    end
    
    if ~isempty(regexpi(shape_struct(i).Name,'[s][h][a][d][o][w]'))
        % Append handles.cad struct (struct or class?)
        handles.cad{end+1}.index = size(data,1);

        %PKZ (20160624) 
        %The number of fields in shape_struct created by
        %LayoverLeopard, effective 20160624, is 9.
        %Number of fields in shape_struct created by RCSTool SaveShapes is 4.
        %The extra 5 fields from LayoverLeopard deal with the stored CAD model.
        %The '(numel(fieldnames(shape_struct)) == 9)' prevents RCSTool from
        %erroring out when handling a non-LayoverLeopard shape_struct.
        
        if (numel(fieldnames(shape_struct)) == 9)
            %convert lat lon to pixel space
            shadowVertices = point_ground_to_slant(shape_struct(i).CADLLAVertices,handles.mitm_hand.Metadata{handles.mitm_hand.Frame});
            shadowVertices = shadowVertices([2 1],:).';   % convert from [row; column] to [column,row]
            shadowVertices = handles.mitm_hand.native2axescoords(shadowVertices); % convert from image coords to axes coords

            % Load CAD model from shape structure - via PlotRangeProfile
            %handles.cad{end}.vertices = shape_struct(i).CADFullImageVertices;
            handles.cad{end}.vertices = shadowVertices;
            handles.cad{end}.faces = shape_struct(i).CADFaces;
            handles.cad{end}.faceAlpha = shape_struct(i).FaceAlpha;
            handles.cad{end}.color = shape_struct(i).Color;
            handles.cad{end}.edgeColor = shape_struct(i).EdgeColor;
            %handles.cad{end}.visibility = 'On';

            %Create the shadow patch
            shadowH = patch('vertices',shadowVertices,'faces',shape_struct(i).CADFaces,'facecolor',shape_struct(i).Color,...
                'edgecolor',shape_struct(i).EdgeColor,'FaceAlpha',shape_struct(i).FaceAlpha);
            set(shadowH, 'Visible', 'off');
            handles.cad{end}.patch = shadowH;
        end % END of 'numel(fieldnames(shape_struct)) == 9' IF
    end
end
%%% End of CAD overlay code by Peter Zweber - BDK 4/20/2016 %%%

set(handles.ShapeTable,'data',data); % Update table data with new shape

guidata(hObject, handles);

% --- Executes on button press in SaveShapes.
function SaveShapes_Callback(hObject, eventdata, handles)
% hObject    handle to SaveShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

[fname, pathstr] = uiputfile( {'*.mat','MATLAB Files (*.mat)'},'Save ROI File',pathstr);
if ~fname, return; end;  % Cancel selected

%grab data from table
data = get(handles.ShapeTable,'data');
if isempty(data)
    disp('No Shapes Available');
    return;
end

shape_struct = struct('Name',{},'Color',{},'Toggle',{},'Vertices',{});
for i=1:size(data,1)
    % We convert all shapes to polygons in lat/long space, since we don't
    % know which image plane/pixel space they will be loaded into.
    shape_struct(i).Name = data{i,1};
    shape_struct(i).Color = data{i,3};
    shape_struct(i).Toggle = data{i,4};
    switch class(handles.rois{i}) % Shape type
        case 'imrect'
            pos = handles.rois{i}.getPosition();
            pos = [pos(1), pos(2); pos(1)+pos(3), pos(2);...
                pos(1)+pos(3), pos(2)+pos(4); pos(1), pos(2)+pos(4)];
        case 'impoly'
            pos = handles.rois{i}.getPosition();
        case 'imellipse'
            pos = handles.rois{i}.getVertices();
    end
    %convert pixels to lat lon space
    for j = 1:size(pos,1)
        pos(j,:) = handles.mitm_hand.axescoords2native(pos(j,:));
    end
    pos = pos(:,[2 1]).'; % convert from [column,row] to [row; column]
    shape_struct(i).Vertices = point_slant_to_ground(pos,...
        handles.mitm_hand.Metadata{handles.mitm_hand.Frame});
end
save(fullfile(pathstr,fname), 'shape_struct');


% --- Executes on button press in dBCheck.
function dBCheck_Callback(hObject, eventdata, handles)
% hObject    handle to dBCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(~isempty(get(handles.ShapeTable,'Data')))
    ComputeRCSData(handles);
end


% --- Executes on button press in SaveKML.
function SaveKML_Callback(hObject, eventdata, handles)
% hObject    handle to SaveKML (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%write KML file that contains RCS results
OverallImageNames = [];

%get last path
if ispref('matlab_sar_toolbox','last_used_directory')
    startpath = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(startpath)||~exist(startpath,'dir')
        startpath = pwd;
    end
else
    startpath = pwd;
end

[fname, path] = uiputfile( {'*.kmz','KMZ Files (*.kmz)'},'Save KMZ File',startpath);
filename = strcat(path,fname);

%grab data from table to get names and use toggle status
data = get(handles.ShapeTable,'data');
if (size(data,1) > 0)
    [pathstr, name] = fileparts(filename);
    filenametemp = sprintf('%s/%s.kml',pathstr,name);
    fid = OpenShapefile(filenametemp,'kml');    
else
    disp('No Shapes Available');
    return;
end

%get directory where we will place image files.  These will get zipped up
%into a KMZ file at the end
pathstr = fileparts(filename);

%compute RCS Data
[SlowData, FastData, SlowXData, FastXData, RangeData, RcsFlags] = ComputeRCSData(handles);  % Added Dr. Zweber RangeData to output list - BDK 4/20/2016

Colors = {'Red','Yellow','Blue','Green','Magenta','White'};    
meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
for i=1:size(data,1)
    Name = data{i,1};
    
    if data{i,4} % Toggle
        %since these shapes may be loaded up in an image with a different
        %azimuth, we will save them all as polygons, sine you can't draw an
        %imrect or imellipse at an arbitrary angle
        switch class(handles.rois{i}) % Shape type
            case 'imrect'
                pos = handles.rois{i}.getPosition();
                pos = [pos(1), pos(2); pos(1)+pos(3), pos(2);...
                    pos(1)+pos(3), pos(2)+pos(4); pos(1), pos(2)+pos(4)];
            case 'impoly'
                pos = handles.rois{i}.getPosition();
            case 'imellipse'
                pos = handles.rois{i}.getVertices();
        end
        %convert pixels to lat lon space
        for j = 1:size(pos,1)
            pos(j,:) = handles.mitm_hand.axescoords2native(pos(j,:));
        end
        llas = point_slant_to_ground(pos.',meta);

        %format shape structure and add to kml file
        Color = Colors{mod(i,6)+1};
        shape.ShapeType = 'POLYGON';
        shape.Lat = llas(1,:);
        shape.Lon = llas(2,:);
        shape.NumVerticies = size(pos,1);
        shape.LineThickness = 2;
        shape.Color = Color;
        shape.FillFlag = 0;
        shape.Count = i*2-1;
        shape.Name = Name;
        AddShape(fid,'kml',shape)
    
        %now we want to add a point at the centroid that contains the image
        %chip (with the shape), some RCS data and the plot of the slow and fast
        %time response curves
        h = PlotRCSData(handles,SlowData,FastData,SlowXData,FastXData,i,RcsFlags);
        
        %save figure as jpg
        PlotJPG = sprintf('%s/%s_Plots.jpg',pathstr,Name);
        f = getframe(h);
        imwrite(f.cdata,PlotJPG,'jpg');        
        close(h); 
        
        %get image region surrounding shape and save as jpg
        MinX = min(pos(:,1));
        MaxX = max(pos(:,1));
        MinY = min(pos(:,2));
        MaxY = max(pos(:,2));
        DeltaX = MaxX-MinX;
        DeltaY = MaxY-MinY;
        CenX = (MaxX+MinX)/2;
        CenY = (MaxY+MinY)/2;
        
        XSize = DeltaX*1.5;
        YSize = DeltaY*1.5;
        
        Size = max([XSize YSize]);
        
        XStart = round(CenX - Size/2);
        YStart = round(CenY - Size/2);
        XStop = round(CenX + Size/2);
        YStop = round(CenY + Size/2);
        
        XStart = max(XStart,1);
        YStart = max(YStart,1);
        XStop = min(XStop,meta.ImageData.NumCols);
        YStop = min(YStop,meta.ImageData.NumRows);

        chip = handles.mitm_hand.readerobj{handles.mitm_hand.Frame}.read_chip...
            ([XStart XStop],[YStart YStop]);
        chip = chip(:,:,1); % Handle multi-channel (polarimetric) case
        h = figure;
        set(h,'Position',[50 50 350 350]);
        
        imshow(densityremap(chip.'),'border','tight'); 
            
        box on;
        %draw shape on image chip
        %compute offset between image on RCSTool screen and image chip 
        OffsetX = XStart;
        OffsetY = YStart;
        ShapeX = pos(:,1) - OffsetX;
        ShapeY = pos(:,2) - OffsetY;        
        
        polyh = impoly(gca,[ShapeX ShapeY]);             
        setColor(polyh,Color);       
        
        f = getframe(h);
        
        %save figure as jpg
        ChipJPG = sprintf('%s/%s_Chip.jpg',pathstr,Name);
        imwrite(f.cdata,ChipJPG,'jpg');        
        close(h);
        
        temp{1} = PlotJPG;
        OverallImageNames = horzcat(OverallImageNames,temp);
        temp{1} = ChipJPG;
        OverallImageNames = horzcat(OverallImageNames,temp);
                        
        Lat = mean(llas(1,:));
        Lon = mean(llas(2,:));        
        
        %Add Placemark
        PlaceName = sprintf('%s_Placemark',Name);
        shape.ShapeType = 'POINT';
        shape.Lat = Lat;
        shape.Lon = Lon;
        shape.Color = Color;
        shape.Count = i*2;
        shape.Name = PlaceName;
        %set description to have formatted html with plots image and data
        PlotShort = sprintf('%s_Plots.jpg',Name);
        ChipShort = sprintf('%s_Chip.jpg',Name);
        shape.Description = GetRCSDescription(meta,Lat,Lon,PlotShort,ChipShort);
        AddShape(fid,'kml',shape);            
    end    
    
    clear ShapeX;
    clear ShapeY;
    clear temp;
end
CloseShapefile(fid,'kml');

%zip kml file and images and rename to kmz file
outfilecell{1} = filenametemp;
OverallImageNames = horzcat(OverallImageNames,outfilecell);
path = fileparts(filename);
tempname = sprintf('%s/temp.zip',path);
zip(tempname,OverallImageNames);
movefile(tempname,filename);
for i=1:length(OverallImageNames)
    delete(OverallImageNames{i});
end

function Description = GetRCSDescription(meta,Lat,Lon,PlotJPG,ChipJPG)

Description = sprintf('<![CDATA[\n');
Description = sprintf('%s<table border="1">\n<tr>\n',Description);
Description = sprintf('%s<td><img src="%s" width="350" height="350" /></td>\n',Description,ChipJPG);
Description = sprintf('%s<td><table border="1" width="444">\n',Description);
Description = sprintf('%s<tr><td><font size="4">IID:</font></td><td><font size="4">%s</font></td></tr>\n',Description,...
                      meta.CollectionInfo.CoreName);
Description = sprintf('%s<tr><td><font size="4">Centroid Lat:</font></td><td><font size="4">%8.5f</font></td></tr>\n',...
                      Description,Lat);
Description = sprintf('%s<tr><td><font size="4">Centroid Lon:</font></td><td><font size="4">%8.5f</font></td></tr>\n',...
                      Description,Lon);
Description = sprintf('%s<tr><td><font size="4">Graze (deg):</font></td><td><font size="4">%5.2f</font></td></tr>\n',Description,...
                      meta.SCPCOA.GrazeAng);
if isfield(meta.SCPCOA,'AzimAng')
    Description = sprintf('%s<tr><td><font size="4">Azimuth (deg):</font></td><td><font size="4">%5.2f</font></td></tr>\n',Description,...
                      meta.SCPCOA.AzimAng);
end
if isfield(meta.CollectionInfo.RadarMode,'ModeID')
      Description = sprintf('%s<tr><td><font size="4">Collect Mode:</font></td><td><font size="4">%s</font></td></tr>\n',Description,...
                    meta.CollectionInfo.RadarMode.ModeID);
end
if isfield(meta.CollectionInfo,'Classification')
    Description = sprintf('%s<tr><td><font size="4">Classification:</font></td><td><font size="4">%s</font></td></tr>\n',Description,...
        meta.CollectionInfo.Classification);
end
Description = sprintf('%s</table></td>\n',Description);
Description = sprintf('%s</table></td>\n',Description);
Description = sprintf('%s</tr></table>\n',Description);
Description = sprintf('%s<table border="1">\n<tr>\n',Description);
Description = sprintf('%s<td><img src="%s" width="800" height="500" /></td>\n',Description,PlotJPG);
Description = sprintf('%s</tr></table>]]>\n',Description);

% --- Executes on button press in SaveRCSData.
function SaveRCSData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveRCSData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get last path
if ispref('matlab_sar_toolbox','last_used_directory')
    startpath = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(startpath)||~exist(startpath,'dir')
        startpath = pwd;
    end
else
    startpath = pwd;
end

%get output filename
[fname, path] = uiputfile( {'*.csv','CSV Files (*.csv)'},'Save CSV File',startpath);

%Compute RCS Data
[SlowData,FastData,SlowXData,FastXData,RangeData,RcsFlags] = ComputeRCSData(handles);

data = get(handles.ShapeTable,'data');
fid = fopen(fullfile(path,fname),'w');

%Write Heading columns with shape name and column titles
for i=1:size(data,1)
    fprintf(fid,'%s,,,,',data{i,1});
end
fprintf(fid,'\n');
MaxSize = 0;
for i=1:size(data,1)
    fprintf(fid,'%s,%s,%s,%s,',RcsFlags.axis_labels.x.slow, ...
            RcsFlags.axis_labels.y,RcsFlags.axis_labels.x.fast, ...
            RcsFlags.axis_labels.y);
    if size(SlowData{i},1) > MaxSize; MaxSize = size(SlowData{i},1); end; 
    if size(FastData{i},1) > MaxSize; MaxSize = size(FastData{i},1); end;
end
fprintf(fid,'\n');

for j=1:MaxSize
    for i=1:size(data,1)
        if size(SlowData{i},1) >= j
            fprintf(fid,'%f,%f,',SlowXData{i}(j),SlowData{i}(j));       
        else
            fprintf(fid,',,');
        end  
        if size(FastData{i},1) >= j
            fprintf(fid,'%f,%f,',FastXData{i}(j),FastData{i}(j));       
        else
            fprintf(fid,',,');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);


function saveShapesNativeCoords(src)

handles = guidata(src);
handles.rois_native = {};
for i = 1:numel(handles.rois)
    handles.rois_native{i}.type = class(handles.rois{i});
    handles.rois_native{i}.color = handles.rois{i}.getColor();
    local_coords = handles.rois{i}.getPosition();
    switch handles.rois_native{i}.type
        case {'imrect', 'imellipse'}
            native_coords = handles.mitm_hand.axescoords2native(local_coords(1:2));
            native_coords(3:4) = local_coords(3:4)*max(handles.mitm_hand.Zoom,1);
        case 'impoly'
            native_coords = zeros(size(local_coords));
            for j = 1:size(local_coords,1)
                native_coords(j,:) = handles.mitm_hand.axescoords2native(local_coords(j,:));
            end
    end
    handles.rois_native{i}.position = native_coords;
end
guidata(src,handles);
    
function restoreShapesLocalCoords(src)

oldzoom = get(zoom(src),'Enable');
oldpan = get(pan(src),'Enable');

handles = guidata(src);
mode = get(get(handles.figure1,'ModeManager'),'CurrentMode');
if verLessThan('matlab','9.5') % MATLAB 2018b on onward
    blocking = ~isempty(mode)&&get(mode,'blocking');
else
    blocking = ~isempty(mode)&&mode.Blocking;
end
for i=1:numel(handles.rois)
    delete(handles.rois{i});
end
handles.rois = {};
for i = 1:numel(handles.rois_native)
    switch handles.rois_native{i}.type
        case {'imrect','imellipse'}
            local_coords = handles.mitm_hand.native2axescoords(...
                handles.rois_native{i}.position(1:2));
            local_coords(3:4) = handles.rois_native{i}.position(3:4)/max(handles.mitm_hand.Zoom,1);
        case 'impoly'
            local_coords = zeros(size(handles.rois_native{i}.position));
            for j = 1:size(local_coords,1)
                local_coords(j,:) = handles.mitm_hand.native2axescoords(...
                    handles.rois_native{i}.position(j,:));
            end
            % If this call is initiated by a zoom-in or zoom-out from a
            % mouse click, then an impoly call will crash here, as will the
            % sets of the zoom and pan states at the end of this function.
            % This is because the zoom state is set to an uniterruptible
            % "on".  The impoly (even when called in a non-interactive
            % way), zoom, and pan functions all attempt to alter the uimode
            % state. Note: imrect and imellipse don't seem to be affected
            % by the blocking state.
            % This got a big messy because in 2018b, Mathworks changed how
            % the stored the blocking state.
            if blocking
                % We temporarily unblock, realizing could be dangerous.
                if verLessThan('matlab','9.5') % MATLAB 2018b on onward
                    set(mode,'blocking',false);
                else
                    mode.Blocking = false;
                end
            end
    end
    handles.rois{i} = feval(handles.rois_native{i}.type, ...
        handles.mitm_hand.AxesHandle, local_coords);
    if blocking  % Return to old state
        if verLessThan('matlab','9.5') % MATLAB 2018b on onward
            set(mode,'blocking',true);
        else
            mode.Blocking = true;
        end
    end
    handles.rois{i}.setColor(handles.rois_native{i}.color);
    handles.rois{i}.addNewPositionCallback(@(pos) ComputeRCSData(guidata(src)));
end
handles.rois_native = {};
guidata(src,handles);

% Restore old state which was changed when we redrew the shapes
try % Doesn't work when initiated from zoom callback.  See comments above.
    zoom(src,oldzoom);
    pan(src,oldpan);
end


function target_az_Callback(hObject, eventdata, handles)
% hObject    handle to target_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_az as text
%        str2double(get(hObject,'String')) returns contents of target_az as a double


% --- Executes during object creation, after setting all properties.
function target_az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in userCalConst.
function userCalConst_Callback(hObject, eventdata, handles)
% hObject    handle to userCalConst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of userCalConst

if get(hObject,'Value')
    set(handles.useCalConstVal,'Enable','on');
    set(handles.db_text,'Enable','on');
else
    set(handles.useCalConstVal,'Enable','off');
    set(handles.db_text,'Enable','off');
end


function useCalConstVal_Callback(hObject, eventdata, handles)
% hObject    handle to useCalConstVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of useCalConstVal as text
%        str2double(get(hObject,'String')) returns contents of useCalConstVal as a double


% --- Executes during object creation, after setting all properties.
function useCalConstVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to useCalConstVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zeropad.
function zeropad_Callback(hObject, eventdata, handles)
% hObject    handle to zeropad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zeropad


function ShapeName_Callback(hObject, eventdata, handles)
% hObject    handle to ShapeName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ShapeName as text
%        str2double(get(hObject,'String')) returns contents of ShapeName as a double


% --- Executes on button press in profile.
function profile_Callback(hObject, eventdata, handles)
% hObject    handle to profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[SlowData, SlowXData, FastData, FastXData, RangeData, RcsFlags]...
    = ComputeRCSData(handles);  %%% Added function call - BDK 4/20/2016
handles.fig_H_RP = PlotRCSRangeProfile(handles,RangeData);  % - BDK 4/21/2016
% Update handles structure
guidata(hObject, handles);   % Added to keep track of RCS plot figure handles - BDK 4/21/2016


function FigHand = PlotRCSRangeProfile(handles,RangeData)
% PlotRangeProfile: Plot RangeProfile (aka Cross Range Sum)
%  Same total result as sum of SlowTime calculation in RCS_Compute

if(~isempty(get(handles.ShapeTable,'Data')))   % Error check to see if data is available in the shape table - BDK 4/21/2016
    if(~isfield(handles,'fig_H_RP'))             % Add field to handles if it doesn't exist - BDK 4/21/2016
        handles.fig_H_RP = [];                     % Initialize if it doesn't exist - BDK 4/21/2016
    end
    
    % Get metadata
    meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
    [cal_sf, cal_bool] = get_cal_factor(handles);

    % Get current Range profile plot figure handles - BDK on 4/21/2016 -- Modified to only use good figure handles 6/8/2016
    FigHand = handles.fig_H_RP(ishandle(handles.fig_H_RP));
    
    % Get shape table data
    data = get(handles.ShapeTable,'data');
    indices = find([data{:,4}]); % Use only the shapes checked "Use"
    
    for i=1:numel(indices)
        h = figure;                       % Create figure - BDK 4/21/2016
        FigHand = [FigHand h];            % Build figure handle vector - BDK 4/21/2016
        set(h,'Name',['Range Profile: ' data{indices(i),1}]); % Added annotation name to Figure - BDK on 4/21/2016
        
        % Plot Range Profile
        RangeYData = meta.Grid.Row.SS*(0:(size(RangeData{i},1)-1));
        plot(RangeData{i},RangeYData);
        
        set(gca,'YDir','reverse');
        title(gca, ['Range Profile: ' data{indices(i),1}]);
        
        % Set ymin, ymax, xmax; xmin calculated below
        ymin = min(RangeYData(:));
        ymax = max(RangeYData(:));
        xmax = max(RangeData{i}(:));
        
        % Set xmin using non -Inf values
        mod_xmin = min(min(RangeData{i}(isfinite(RangeData{i}))));
        
        % Set axes range
        axis([mod_xmin-1 xmax+1 ymin ymax]);
        
        % Assign axes labels
        grid on;
        measure_type_strings = get(handles.measuretype,'String');
        xlabel(gca,measure_label(...
            measure_type_strings{get(handles.measuretype,'Value')}, ...
            get(handles.dBCheck, 'Value'), cal_bool));
        ylabel(gca,'Range from top of ROI [m]');
        
        % Modify and plot CAD model if it exists
        if ~isempty(handles.cad)
            for jj=1:numel(handles.cad)
                if handles.cad{jj}.index == indices(i)
                    
                    %PKZ (20160624)
                    %If a shadow or layover annotation has an accompanying
                    %CAD model, numel(fieldnames(handles.cad{jj})) == 7),
                    %if it doesn't have an accompanying CAD model,
                    %do not attempt to overlay the non-existent model on the
                    %figure.  duh!
                    if (numel(fieldnames(handles.cad{jj})) ~= 7)
                        continue;
                    end
                    cad = handles.cad{jj};
                    layVertRP = cad.vertices;
                    
                    % Find ymin and ymax of CAD model vertices
                    v_ymin = min(layVertRP(:,2));
                    v_ymax = max(layVertRP(:,2));
                    v_ylen = v_ymax - v_ymin;
                    
                    % Recenter vertices y-coordinates to y=0
                    layVertRP(:,2) = layVertRP(:,2) - v_ymin;
                    
                    % Rescale vertices
                    rpAxesY = get(gca,'YLim');
                    ysf = (rpAxesY(2)-rpAxesY(1))/v_ylen;
                    layVertRP(:,2) = layVertRP(:,2)*ysf;
                    
                    % Recenter vertices x-coordinates to x=0
                    v_xmin = min(layVertRP(:,1));
                    v_xmax = max(layVertRP(:,1));
                    v_xlen = v_xmax - v_xmin;
                    layVertRP(:,1) = layVertRP(:,1) - v_xmin;
                    
                    % Rescale x-dim vertices
                    rpAxesX = [mod_xmin xmax];
                    xsf = (rpAxesX(2)-rpAxesX(1))/v_xlen;
                    layVertRP(:,1) = layVertRP(:,1)*xsf;
                    
                    % Shift so lowest x value of CAD model to plot left
                    layVertRP(:,1) = layVertRP(:,1) + rpAxesX(1);
                    
                    % Draw the layover patch on the range profile
                    handles.layoverRangeProfileH2D = patch('Parent',gca,...
                        'vertices',layVertRP,'faces',cad.faces,'facecolor',cad.color,...
                        'edgecolor',cad.edgeColor,'FaceAlpha',cad.faceAlpha,...
                        'Visible','on');
                    
                    % Replot RangeData on top of CAD model
                    hold on;
                    plot(RangeData{i},RangeYData,'k','LineWidth',1.3);
                    
                    % Rotate range profile and CAD model if available so range is along the x-axis and
                    % amplitude is along the y-axis (PKZ: Removed comment) - BDK 4/21/2016
                    % view(90,-90);  % Commented out plot rotation so range profile is oriented like the SAR image. - BDK 4/21/2016
                    % axis xy;       % Commented out plot rotation so range profile is oriented like the SAR image. - BDK 4/21/2016
                    % ylabel(' Range ---> (m)');   % Commented out plot rotation so range profile is oriented like the SAR image. - BDK 4/21/2016
                    % End of range profile rotation  - BDK 4/21/2016
                    
                end
            end
        end
        % end  % Removed if-statement so plots are always generated - BDK 6/23/2016
    end
end

% --- Executes on selection change in slowtimeunits.
function slowtimeunits_Callback(hObject, eventdata, handles)
% hObject    handle to slowtimeunits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns slowtimeunits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from slowtimeunits

combo_string = get(hObject,'String');
if strcmpi(combo_string{get(hObject,'Value')},'Target Relative')
    set(handles.target_az,'Enable','on');
    set(handles.targetaz_text,'Enable','on');
else
    set(handles.target_az,'Enable','off');
    set(handles.targetaz_text,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function slowtimeunits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slowtimeunits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select All Button - BDK 3/18/2016 %
data = get(handles.ShapeTable,'Data');   % Get all of the shape table data
for r = 1:size(data,1)                   % Select all rows in shape table - BDK 3/18/2016
  data{r,4} = true;
end
set(handles.ShapeTable, 'Data', data);     % Change values in shape table in GUI


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deselect All Button  - BDK 3/18/2016 %
data = get(handles.ShapeTable,'Data');  % Get all of the shape table data
for r = 1:size(data,1)                   % Select all rows in shape table - BDK 3/18/2016
  data{r,4} = false;
end
set(handles.ShapeTable, 'Data', data);    % Change values in shape table in GUI


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Jump/Zoom Button -  BDK 3/24/2016
% Get ROI data from table.
data = get(handles.ShapeTable,'data');
selectedroi = find([data{:,4}]);   % Determine which annotation is selected - BDK 3/24/2016
if isempty(data)
    msgbox('Jump/Zoom - No Shapes Available','RCSTool');  % Modified so text is displayed in a pop-up window to alert the user. - BDK 3/24/2016
    return;
end

if isempty(selectedroi)
    msgbox('Jump/Zoom - No Shape Selected','RCSTool');   % Modified so text is displayed in a pop-up window to alert the user. - BDK 3/24/2016
    return
end

if(length(handles.selectedroi)> 1)   % Check to ensure only one annotation is selected - BDK 3/24/2016
   msgbox('Jump/Zoom - Too Many Shapes Selected','RCSTool');  % Alert user when more than one annotation is selected - BDK 3/24/2016
   return                             % BDK 3/24/2016
end                                   % BDK 3/24/2016

switch class(handles.rois{selectedroi}) % Shape type
    case 'imrect'
        pos = handles.rois{selectedroi}.getPosition();
        pos = [pos(1), pos(2); pos(1)+pos(3), pos(2);...
            pos(1)+pos(3), pos(2)+pos(4); pos(1), pos(2)+pos(4)];
    case 'impoly'
        pos = handles.rois{selectedroi}.getPosition();
    case 'imellipse'
        pos = handles.rois{selectedroi}.getVertices();
end
handles.mitm_hand.setView('Zoom', 1, 'CenterPos', ...
    handles.mitm_hand.axescoords2native(mean(pos)));


% --- Executes on button press in BatchRCS.
function BatchRCS_Callback(hObject, eventdata, handles)
% hObject    handle to BatchRCS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get input directory for files
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end
inputdir = uigetdir(pathstr,'Select Input Folder');

%have user specify output csv filename (also will write an html file, maybe
%kml)
[outfile,outdir] = uiputfile('*.csv','Save Batch Output Data',inputdir);
outfile = [outdir outfile];

%get reference filename and shape information (name, lat/lon coordinates)
data = get(handles.ShapeTable,'data');
indices = find([data{:,4}]); % Use only the shapes checked "Use"

for i=1:numel(indices)
    % Getting shape positions
    switch data{indices(i),2} % Shape type
        case 'Rectangle'
            pos = handles.rois{indices(i)}.getPosition();
            pos = [pos(1), pos(2); pos(1)+pos(3), pos(2);...
                pos(1)+pos(3), pos(2)+pos(4); pos(1), pos(2)+pos(4)];
        case 'Polygon'
            pos = handles.rois{indices(i)}.getPosition();
        case 'Ellipse'
            pos = handles.rois{indices(i)}.getVertices();
    end
    for j = 1:size(pos,1)
        pos(j,:) = handles.mitm_hand.axescoords2native(pos(j,:));
    end
    
    Shapes(i).Name = data{indices(i),1};
    Shapes(i).Pos = pos;
end

RefFile = handles.mitm_hand.Filename{1};

%run Batch RCS
RunBatchRCS(RefFile,Shapes,inputdir,outfile);


% --- Executes on selection change in measuretype.
function measuretype_Callback(hObject, eventdata, handles)
% hObject    handle to measuretype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(~isempty(get(handles.ShapeTable,'Data')))
    ComputeRCSData(handles);
end


% --- Executes during object creation, after setting all properties.
function measuretype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to measuretype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
