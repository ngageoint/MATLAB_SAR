% Load AOI from file
function ApToolLoadImage(filenames,hObject,handles,varargin)

if isfield(handles,'H')
    delete(handles.H);
end

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
    if(iscell(filenames)), filenames=sort(filenames); end
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

[readerobj_indices,reader_obj] = group_reader_objs_by_pol(reader_obj); % Consolidate reader objects into polarimetric sets

%get metadata
if (iscell(reader_obj))
    meta=reader_obj{1}.get_meta(); % Assume first frame for getting datasize
    % Adjust later if frame changes
else
    meta=reader_obj.get_meta();
end
nx = double(meta.ImageData.NumCols);
ny = double(meta.ImageData.NumRows);

%set filename or corename as the title of the GUI
if ~iscell(filenames)
    temp = filenames;
    clear filenames;
    filenames{1} = temp;
end
[~,name,ext] = fileparts(filenames{1});
fname = strcat(name,ext);
set(handles.figure1,'Name',['ApertureTool, Filename: ', fname]);

% Select AOI and frame number if necessary
if isempty(varargin) || isempty(varargin{1}) || varargin{1}(1) == 0
    %determine if we'll just load the image or if we'll load an AOI
    %selection MITM_Viewer    
    if (nx <= handles.MaxLoadSize && ny <= handles.MaxLoadSize)
        %load image directly into tool
        aoi = [1 1 nx ny];
        FrameNum = 1;
    else
        %launch a MITM viewer with an AOI selection box
         ButtonName = questdlg(sprintf('%s %i x %i\n%s %i x %i\n%s',...
            'Suggested AOI dimensions are less than',handles.MaxLoadSize,handles.MaxLoadSize,...
            'Current AOI dimensions are',nx,ny,...
            'Do you want to redefine the AOI?'), ...
            'AptertureTool AOI', ...
            'Yes', 'No', 'Yes');
        if isempty(ButtonName), return; end %window was closed CANCEL
        if strcmpi(ButtonName,'No')
            aoi = [1 1 nx ny]; %use AOI anyway
            FrameNum = 1;
        else
            [aoi, FrameNum] = mitm_viewer(fullfilename,'mode','aoi','closeAfterSelect',true);
            FrameNum = find(cellfun(@(x) isequal(FrameNum,sf_ind(x)), readerobj_indices));
        end
        if isempty(aoi), return; end % Cancel was selected
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

set(handles.nx,'String',nx);
set(handles.ny,'String',ny);

%Draw metaicon
MetaIcon(meta,'handle',handles.metaicon);

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
    set(handles.UniformWeighting,'Value',true);
    set(handles.UniformWeighting,'Enable','off');
else
    set(handles.UniformWeighting,'Value',false);
    set(handles.UniformWeighting,'Enable','on');
end
% Set default values for required frequency support deskew values
handles.isvalid_row_deltakcoapoly = isfield(meta,'Grid') && ...
    isfield(meta.Grid,'Row') && isfield(meta.Grid.Row,'DeltaKCOAPoly') && ...
    ~isempty(meta.Grid.Row.DeltaKCOAPoly) && ...
    all(isfinite(meta.Grid.Row.DeltaKCOAPoly(:)));
if handles.isvalid_row_deltakcoapoly && any(meta.Grid.Row.DeltaKCOAPoly(:)~=0)
    set(handles.DeskewFast,'Enable','on');
else
    set(handles.DeskewFast,'Enable','off');
    meta.Grid.Row.DeltaKCOAPoly = 0; % Give default value (no shift)
end
handles.isvalid_col_deltakcoapoly = isfield(meta,'Grid') && ...
    isfield(meta.Grid,'Col') && isfield(meta.Grid.Col,'DeltaKCOAPoly') && ...
    ~isempty(meta.Grid.Col.DeltaKCOAPoly) && ...
    all(isfinite(meta.Grid.Col.DeltaKCOAPoly(:)));
if handles.isvalid_col_deltakcoapoly && any(meta.Grid.Col.DeltaKCOAPoly(:)~=0)
    set(handles.DeskewSlow,'Enable','on');
else
    set(handles.DeskewSlow,'Enable','off');
    meta.Grid.Col.DeltaKCOAPoly = 0; % Give default value (no shift)
end
% Is it possible to deskew frequency support in both dimensions at once?
if isscalar(meta.Grid.Col.DeltaKCOAPoly) && ...
        isscalar(meta.Grid.Row.DeltaKCOAPoly)
    set(handles.DeskewSlow,'Value',handles.isvalid_col_deltakcoapoly); % Center both dimensions
    set(handles.DeskewFast,'Value',handles.isvalid_row_deltakcoapoly);
else % Make user select dimension to deskew
    set(handles.DeskewSlow,'Value',handles.isvalid_col_deltakcoapoly && ...
        all(meta.Grid.Col.DeltaKCOAPoly(:) == 0));
    set(handles.DeskewFast,'Value',handles.isvalid_row_deltakcoapoly && ...
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
        set(handles.image_rep_combo,'value',1);
    else
        set(handles.image_rep_combo,'value',...
            find(strcmpi('Pauli',get(handles.image_rep_combo,'string'))));
    end
    handles.polar_PHD_idx = 1;
    if numel(handles.TxRcvPolarizationProc) == 1
        set(handles.image_rep_combo,'value',1,'visible','off');
        set(handles.PHDPolCombo,'visible','off');
    else
        set(handles.image_rep_combo,'value',2,'visible','on');
        set(handles.PHDPolCombo,'visible','on');
    end
    %set strings in PHDPolCombo
    set(handles.PHDPolCombo,'String',meta.ImageFormation.TxRcvPolarizationProc);
else
    set(handles.image_rep_combo,'value',1,'visible','off');
    set(handles.PHDPolCombo,'visible','off');
end


% Read data
handles.complex_data = double(reader_obj{FrameNum}.read_chip(...
    [aoi(1) aoi(1)+aoi(3)-1],[aoi(2) aoi(2)+aoi(4)-1]));
reader_obj{FrameNum}.close();
% Condition data so that imcontrast control works better
handles.complex_data = 10*handles.complex_data./mean(handles.complex_data(:));
handles.aoi = aoi; % store aoi for later
handles.segment = FrameNum;

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

%set initial SS, Ap Time and Ap Bandwidth

handles.meta = meta; % store metadata for later

ApToolprocess_phd(hObject,handles); % Apply PHD options
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
api_handle.addNewPositionCallback(@(x) ApToolnewpos(x,hObject));

% Now enable relevant controls
set(handles.Play,'enable','on');
set(handles.Last,'enable','on');
set(handles.First,'enable','on');
set(handles.Next,'enable','off'); % Full aperture is not a frame, so
set(handles.Prev,'enable','off'); % no next or previous
set(handles.Recording,'Enable','on');
set(handles.SaveImage,'Enable','on');
set(handles.ManualAdjust,'enable','on');

% Update handles structure
guidata(hObject, handles);

% Update resolution fields
ApToolnewpos(getPosition(handles.H),hObject);
