classdef hg_mitm_viewer < hgsetget
%HG_MITM_VIEWER MATLAB graphics object for More-Image-Than-Memory complex image viewer
%   H = HG_MITM_VIEWER (PARENT_UIPANEL) displays a file (or files) of
%   complex data in a MATLAB IMAGE/IMAGESC/IMSHOW style viewer, even if the
%   file is too large to fit into memory.  Any file type handled by
%   OPEN_READER can be used.
%
%   The purpose of this function is to allow a user to easily view and
%   navigate complex imagery too large to fit into memory.  It is intended
%   to represent data pixel-for-pixel as truly as possible to the
%   underlying data without any distortion.  There is no projection,
%   filtering, or processing of any sort outside of the user-defined
%   decimation and remaps.
%
%   A data cursor tool is available, just as in a standard MATLAB figure.
%   However, this tool provides the coordinates for each pixel in the full
%   complex image that is stored in the file, not just the portion of the
%   image shown on the screen.  It also shows the value of the pixel at
%   that point, both in complex and magnitude/angle format.
%
%   By default, the figure initially opens with an overview of the whole
%   data file, but the figure supports interactive zooming and panning.
%   Only as much resolution and extent as are needed for the current figure
%   size and zoom level are read from file.  This allows for viewing of
%   extremely large datasets even on computers with limited memory.
%   Zooming can be done with the magnifying glass tools in the toolbar or
%   the scale factor dropbox in the toolbar.  Single click zoom and window
%   zoom should both work.  Panning is done with the hand tool just as in a
%   standard MATLAB figure.  New data will be read from the file after each
%   zoom and pan.
%   
%   Several decimation types and remaps are available.  The dropboxes in
%   the toolbar allow a user to select these.  The dropboxes are editable
%   if the user wants to use his own decimation or remap function. Just
%   type the name of a function that can be found on the MATLAB path in the
%   drop box.  Decimation function operate along columns to produce a
%   single row-- like the MATLAB max and mean function. (Type in "min" to
%   the decimation dropbox to see how an external function can be used to
%   decimated.)  Remap functions take a single image as in input and output
%   a single image as an output.
%
%   This class has a number of properties that can be read and set:
%
%       Property name     Description
%       AxesHandle        MATLAB handle for axes used for displaying data.
%                         Read only.
%       CenterPos         Center of current view in pixel coordinates from
%                         the full image.
%       ComplexData       Raw data (not remapped) in current view (with
%                         decimation).  Read only.
%       DataTransformFcn  Callback that transforms the underlying complex
%                         data before remap and display.
%       Decimation        Type of decimation used for reading and
%                         displaying data.
%       Filename          Full path and filename(s) of currently opened
%                         file(s).  Cell array of strings (in case multiple
%                         files are open.  Read only.
%       Frame             Number of frame currently being viewed for
%                         multi-image datasets.
%       FrameFileIndices  Index to which of the input file(s) passed to
%                         openFile correspond to which frame in the
%                         multiframe view.  Can be multiple.  Read only.
%       FrameSubFileIndices Index to which subimage(s) (from a multi-image
%                           file) correspond to which frame in the
%                           multiframe view.  Can be multiple.  Read only.
%       Metadata          Structure with the SICD description of the data.
%                         Read only.
%       Parent            Uipanel component in which to place viewer.  Read
%                         only (although settable on initialization.)
%       Position          Position of view in full file in vector form:
%                            [left top width height]
%                         where view is assumed to be energy-from-top,
%                         view-from-above.  Read only.  Position can be
%                         changed though by setting CenterPos and Zoom.
%       PreChangeViewFcn  Callback that occurs before any new data is
%                         loaded into view
%       PostChangeViewFcn Callback that occurs after any new data is loaded
%                         into view.
%       Remap             Remap function.  Default is 'densityremap'.
%       Zoom              Zoom level.  2 means one pixel on the screen
%                         represents a 2x2 pixel area in the file.
%
%   There are also a number of public methods that can be called:
%       Method name       Description
%       axescoords2native Convert coordinates in the MATLAB axes to
%                         coordinates in the full image space.
%       close             Close currently open file.
%       geojump           Place data cursor at a given lat/long/altitude.
%       native2axescoords Inverse of axescoords2native.
%       openFile          Open one (or more) files for viewing.
%       setView           Change current view of currently loaded complex
%                         data.  Allows for changing multiple view
%                         properties simultaneously.
%
% WARNING:
% This class only works for MATLAB version 7.7 (2008b) and above.
% The object-oriented features of MATLAB were not introduced until version
% 7.6 (2008a).  It also uses javaObjectEDT (which is much faster than the
% older awtcreate/awtinvoke) which was introduced in 7.7 (2008b).
%
% Note that this function depends heavily on using Java components within
% MATLAB.  This is unsupported and undocumented feature of MATLAB.  This
% may not work in a future version.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Properties
% Constants
properties (Constant, GetAccess = private, Hidden)
    TOOLBAR_HEIGHT = 25; % Default height in pixels for MATLAB toolbars
    DEFAULT_COLORMAP = gray(256); % 8-bit grayscale
    DEFAULT_REMAP = 'densityremap';
    DEFAULT_POL_DECOMP = 'pol_decomp.Pauli';
end

% Read-only
properties (GetAccess = public, SetAccess = private)
    AxesHandle; % MATLAB handle for axes used for displaying data
    Filename;
    % Index to which of the input file(s) passed to openFile correspond to
    % which frame in the multiframe view.  Can be multiple.
    FrameFileIndices;
    % Index to which subimage(s) (from a multi-image file) correspond to
    % which frame in the multiframe view.  Can be multiple.
    FrameSubFileIndices;
    Metadata; % SICD metadata structure of current dataset
    Parent; % Parent uipanel
end

properties (Dependent = true, GetAccess = public)
    % Read-only
    ComplexData; % Raw data (not remapped) in current view (with decimation)
    Position; % Position of view in full file [left top width height]
    
    % Public read/write
    CenterPos; % Center of current view in pixel coordinates from the full image
    Decimation; % Decimation function
    Frame; % Frame we are currently viewing
    Remap; % Remap used to go from complex data to pixel value
    Zoom; % Zoom (or decimation) level
end

% No special get/set necessary, so just make public access
properties (GetAccess = public, SetAccess = public)
    PreChangeViewFcn = @deal; % Callback that occurs before any new data is loaded into view
    PostChangeViewFcn = @deal; % Callback that occurs after any new data is loaded into view
    DataTransformFcn = []; % Callback that transforms the underlying complex data before remap and display
end

% Want to allow access for advanced use, but not for general use, so we hide
properties (GetAccess = public, SetAccess = private, Hidden)
    main_toolbar; % Java handle to main toolbar
    movie_toolbar; % Java handle to movie toolbar
    readerobj; % Reader object(s) from open_reader.m for current data
end

% Properties only used internally
properties (GetAccess = private, SetAccess = private, Hidden)
    % GUI components
    figure_handle; % MATLAB handle for parent figure
    image_handle; % MATLAB handle for image object used for displaying data
    main_toolbar_buttons; % Java handles to buttons on main toolbar
    main_toolbar_hg; % MATLAB handle graphics container for main toolbar
    movie_toolbar_hg; % MATLAB handle graphics container for movie toolbar
    
    % Toolbar component Java handles
    decimateComboBox;
    frameComboBox;
    remapComboBox;
    zoomComboBox;

    % Image display properties
    climlist; % List of colormap limits for each frame (for floating point remaps)
    cmaplist; % List of colormaps for each frame
    remaplist; % List of remaps for each frame
    decimation_function = 'none'; % Current decimation function (one for all frames)
    frames_loaded; % Array of which frames currently have the displayed area in memory

    % File/data properties
    complex_image = {}; % Actually complex data for the currently displayed pixels

    % Navigation properties
    current_frame = 1; % Number of frame we are currently viewing
    movie_timer; % Timer object that supports animating multiple frames
    origin = [1 1]; % Pixel coordinate of upper-left corner of displayed area
    panzoom = false; % Flag to prevent multiple concurrent pan/zooms
    pan_timer; % Timer object that supports panning
    scale = 1; % Zoom (or decimation) level
    axes_size_prezoom = [];
end

% Dependent properties.  Used only for each of notation
properties (Dependent = true, GetAccess = private, Hidden)
    axes_size; % Size in pixels of the axes being used
    data_size; % Size in pixels of entire file.
    image_size; % Size of image currently displayed (slightly different
                % than the axes size for reasons we won't go into here)
    number_of_frames; % Total number of frames currently available
    total_toolbar_height; % Height in pixels user by toolbar(s)
end

methods
    %% Constructor
    function obj = hg_mitm_viewer(parent_uipanel)
        % Check validity of input argument
        if ~ishghandle(parent_uipanel)||~strcmpi(get(parent_uipanel,'type'),'uipanel')
            error('HG_MITM_VIEWER:InvalidParent','Parent handle must be a valid uipanel.');
        end
        
        % Setup figure layout
        set(parent_uipanel, 'BorderType', 'none');
        obj.Parent = parent_uipanel;
        set(obj.Parent,'DeleteFcn',@(src,event) obj.delete());
        obj.figure_handle = ancestor(parent_uipanel,'figure');
        set(obj.figure_handle,'Colormap',obj.DEFAULT_COLORMAP);
        % MJToolBar is the MATLAB wrapper for the Java JToolBar.  It is useful
        % because it adds the overflow icon when the toolbar becomes too small to
        % hold all of its components.
        obj.movie_toolbar = javaObjectEDT('com.mathworks.mwswing.MJToolBar');
        obj.movie_toolbar.setFloatable(false); % Gets rid of goofy grab "dots" on left of toolbar
        [movie_toolbar_java,obj.movie_toolbar_hg]=javacomponent(obj.movie_toolbar,[],obj.Parent);
        obj.movie_toolbar.setVisible(false); % Not visible until we have multiple frames
        obj.add_movie_toolbar_buttons();
        % We add the movie toolbar first, so that if it is for some reason
        % made visible, it still falls behind the main toolbar.
        obj.main_toolbar = javaObjectEDT('com.mathworks.mwswing.MJToolBar');
        obj.main_toolbar.setFloatable(false); % Gets rid of goofy grab "dots" on left of toolbar
        [main_toolbar_java,obj.main_toolbar_hg]=javacomponent(obj.main_toolbar,[],obj.Parent);
        obj.add_main_toolbar_buttons();
        obj.movie_timer = timer('ExecutionMode', 'fixedRate', ...
            'TimerFcn', @(hco, user) obj.timerTickFcn(), ...
            'StopFcn', @(hco, user) obj.timerStopFcn(), ...
            'BusyMode', 'drop', 'TasksToExecute', inf, 'Period', 0.5);
        obj.AxesHandle = axes('Parent',obj.Parent,'Visible','off','Position',[0 0 1 1]);
        obj.image_handle = image(0,'Parent',obj.AxesHandle,'CDataMapping','scaled'); % Create empty MATLAB image object for later use
        set(obj.AxesHandle,'Visible','off','DataAspectRatio',[1 1 1],'CLim',[0 1]); % image() changes these

        % Setup zoom scale and zoom scale dropbox
        obj.zoomComboBox=javaObjectEDT('javax.swing.JComboBox');
        obj.zoomComboBox.setToolTipText('Zoom Level');
        obj.zoomComboBox.setMinimumSize(javaObjectEDT('java.awt.Dimension',60,16));
        obj.zoomComboBox.setMaximumSize(javaObjectEDT('java.awt.Dimension',60,1000));
        obj.zoomComboBox.setEditable(true);
        set(handle(obj.zoomComboBox,'callbackproperties'),...
            'actionPerformedCallback',@(src,event) updateZoomComboBox(obj,event));

        % Setup decimation types and dropbox
        if(~isempty(ver('images'))) % Image Processing toolbox require for decimation feature in readers
            obj.decimateComboBox=javaObjectEDT('javax.swing.JComboBox',{'Decimate','max','mean'});
            obj.decimateComboBox.setToolTipText('Decimation Type');
            obj.decimateComboBox.setMinimumSize(javaObjectEDT('java.awt.Dimension',100,16));
            obj.decimateComboBox.setMaximumSize(javaObjectEDT('java.awt.Dimension',100,1000));
            obj.decimateComboBox.setEditable(true);
            set(handle(obj.decimateComboBox,'callbackproperties'),...
                'actionPerformedCallback',@(src,event) updateDecimateComboBox(obj,event));
        end

        % Setup remap types and dropbox
        obj.remapComboBox=javaObjectEDT('javax.swing.JComboBox');
        obj.populateRemapComboBox();
        obj.remapComboBox.setToolTipText('Remap Type');
        obj.remapComboBox.setMinimumSize(javaObjectEDT('java.awt.Dimension',125,16));
        obj.remapComboBox.setMaximumSize(javaObjectEDT('java.awt.Dimension',125,1000));
        obj.remapComboBox.setEditable(true);
        obj.remapComboBox.setSelectedItem(obj.DEFAULT_REMAP);
        set(handle(obj.remapComboBox,'callbackproperties'),...
            'ActionPerformedCallback',@(src,event) updateRemapComboBox(obj,event));

        % Setup frame selector
        obj.frameComboBox=javaObjectEDT('javax.swing.JComboBox');
        obj.frameComboBox.setToolTipText('Frame Number');
        obj.frameComboBox.setMinimumSize(javaObjectEDT('java.awt.Dimension',50,16));
        obj.frameComboBox.setMaximumSize(javaObjectEDT('java.awt.Dimension',50,1000));
        set(handle(obj.frameComboBox,'callbackproperties'),...
            'ActionPerformedCallback',@(src,event) updateFrameComboBox(obj,event));
        
        % Draw Java combobox components on toolbars
        obj.main_toolbar.addSeparator;
        obj.main_toolbar.add(obj.zoomComboBox);
        if(~isempty(ver('images'))), obj.main_toolbar.add(obj.decimateComboBox); end;
        obj.main_toolbar.add(obj.remapComboBox);
        obj.main_toolbar.repaint;
        obj.main_toolbar.revalidate;
        obj.movie_toolbar.addSeparator;
        obj.movie_toolbar.add(obj.frameComboBox);
        obj.movie_toolbar.repaint;
        obj.movie_toolbar.revalidate;

        % Setup figure toolbar callbacks
        set(datacursormode(obj.figure_handle),'UpdateFcn',@(src, event) obj.complexdatacursorfun(event));
        set(zoom(obj.figure_handle),'ButtonDownFilter',@(src, event) obj.noZoomOut);
        set(zoom(obj.figure_handle),'ActionPostCallback',@(src, event) obj.handleZoomClick(src));
        set(pan(obj.figure_handle),'ActionPreCallback',@(src, event) obj.prePanAction);
        set(pan(obj.figure_handle),'ActionPostCallback',@(src, event) obj.handlePanClick);
        obj.pan_timer=timer('TimerFcn',@(src,event) obj.handlePanClickTimer,'StartDelay',.1);

        % Layout manager
        obj.mitmResizeFcn(); % Call once in case no resize events
        set(obj.Parent,'ResizeFcn',@(src,event) mitmResizeFcn(obj));
    end % Constructor
    
    %% Destructor
    function delete(obj)
        close(obj);
        delete(obj.pan_timer);
        delete(obj.movie_timer);
    end % Destructor

    %% Open file
    function openFile(obj,filenames,suppress_load)
        if nargin<3 % Allow for opening file, but not loading any pixel data
            suppress_load = false;
        end
        
        obj.close(); % Close any previously open files
        
        if ischar(filenames) || isinmem(filenames) % Input can be single dataset or cell array of datasets
            filenames = {filenames}; % Just treat both types as cell array
        end
        obj.Filename = filenames;
        % Open readers for all files, and build cell array of all readers
        sf_ind = []; fn_ind = []; % Indices to filenames and subimages for each reader object
        for i = 1:length(filenames)
            newreader = open_reader(filenames{i});
            fn_ind = [fn_ind i*ones(size(newreader))];
            sf_ind = [sf_ind 1:numel(newreader)];
            obj.readerobj = [obj.readerobj newreader];
        end
        clear newreader;
        % Some images are spread across multiple files
        [readerobj_indices,obj.readerobj]=group_reader_objs_by_pol(obj.readerobj); % Consolidate reader objects into polarimetric sets
        obj.FrameFileIndices = {}; obj.FrameSubFileIndices = {};
        for i = 1:numel(readerobj_indices)
            obj.FrameFileIndices{i} = fn_ind(readerobj_indices{i});
            obj.FrameSubFileIndices{i} = sf_ind(readerobj_indices{i});
        end
        obj.populateFrameComboBox();
        % Extract metadata for each image
        obj.Metadata=cell(size(obj.readerobj));
        for i=1:obj.number_of_frames
            obj.Metadata{i}=obj.readerobj{i}.get_meta();
            if((obj.Metadata{i}.ImageData.NumRows~=obj.Metadata{1}.ImageData.NumRows)||...
                    (obj.Metadata{i}.ImageData.NumCols~=obj.Metadata{1}.ImageData.NumCols))
                warning('mitm_viewer:nonmatching_sizes','Sizes of all images passed to mitm_viewer are not the same.  Some features of viewer might not work.');
                continue;
            end
        end
        obj.populateZoomComboBox();
        
        % Setup multi-frame memory of display parameters
        obj.frames_loaded = false(obj.number_of_frames,1);
        obj.climlist = cell(obj.number_of_frames,1);
        obj.cmaplist = cell(obj.number_of_frames,1);
        obj.remaplist = cell(obj.number_of_frames,1);
        for i = 1:obj.number_of_frames
            obj.cmaplist{i} = obj.DEFAULT_COLORMAP;
            obj.remaplist{i} = get(obj.remapComboBox,'SelectedItem');
        end
        
        % Display full image extent
        obj.mitmResizeFcn(); % Draw/remove movie toolbar if necessary
        if ~suppress_load
            max_img_size = getpixelposition(obj.Parent);
            max_img_size = max_img_size([3 4])-[0 obj.total_toolbar_height]; % Minus toolbars
            max_img_size = max(max_img_size,1); % Parent could be too small
            initial_zoom = max(ceil(obj.data_size./max_img_size)); % Initial scale is lowest integer decimation that will fit
            obj.setView('Zoom', initial_zoom, 'CenterPos', [1 1]);
        end
    end
    
    %% Close file
    function close(obj)
        % Clear timers
        if ~isempty(obj.pan_timer)||~isempty(obj.movie_timer)
            while(strcmp(get(obj.pan_timer,'Running'),'on'));
                stop(obj.pan_timer);
            end
            while(strcmp(get(obj.movie_timer,'Running'),'on'));
                stop(obj.movie_timer);
            end
        end
        % Close readers if open
        if ~isempty(obj.readerobj)
            for i=1:obj.number_of_frames
                obj.readerobj{i}.close();
                obj.readerobj{i} = [];
            end
        end
        % Reset private properties
        obj.readerobj = {};
        obj.complex_image = {};
        set(obj.image_handle,'CData',0);
        obj.current_frame = 1;
        obj.decimation_function = 'none';
        obj.remaplist = [];
        obj.scale = 1;
        % Clear all public properties as well
        % For now we only do a few.  We probably should to all.
        obj.Filename = [];
        obj.FrameFileIndices = [];
        obj.FrameSubFileIndices = [];
        obj.Metadata = [];
    end
    
    %% Get/Set public methods
    function coords = get.CenterPos(obj)
        coords = obj.axescoords2native([...
            mean(get(obj.AxesHandle,'Xlim'))...
            mean(get(obj.AxesHandle,'Ylim'))]...
            -((obj.axes_size-obj.image_size)/2)); % Due to the fact that the axes is slightly larger than image
    end
    function set.CenterPos(obj, val)
        obj.setView('CenterPos', val);
    end
    
    function ComplexData = get.ComplexData(obj)
        if numel(obj.complex_image)>=obj.current_frame % Might not have opened any data yet
            ComplexData = obj.complex_image{obj.current_frame};
        else
            ComplexData = [];
        end
    end
    
    function set.DataTransformFcn(obj,val)
        obj.DataTransformFcn = val;
        obj.applyNewRemapOnly(obj.Remap);
    end

    function Decimation = get.Decimation(obj)
        Decimation=obj.decimation_function;
    end
    function set.Decimation(obj, val)
        obj.setView('Decimation', val);
    end

    function Frame = get.Frame(obj)
        Frame = obj.current_frame;
    end
    function set.Frame(obj, val)
        obj.setView('Frame', val);
    end
    
    function Position = get.Position(obj)
        if(obj.scale<1)
            screenxlim=get(obj.AxesHandle,'Xlim');
            screenylim=get(obj.AxesHandle,'Ylim');
            Position = [obj.origin + [screenxlim(1) screenylim(1)] - 1 ... % obj.origin is the origin of the 1:1 area, not the current zoom
                obj.image_size*obj.scale];
        else
            Position = [obj.origin obj.image_size*obj.scale];
        end
        Position(3:4) = min(Position(3:4),obj.data_size);
    end
    
    function Remap = get.Remap(obj)
        if obj.current_frame<=numel(obj.remaplist)
            Remap = obj.remaplist{obj.current_frame};
        else
            Remap = obj.DEFAULT_REMAP;
        end
    end
    function set.Remap(obj, val)
        obj.setView('Remap', val);
    end
    
    function Zoom = get.Zoom(obj)
        Zoom = obj.scale;
    end
    function set.Zoom(obj, val)
        obj.setView('Zoom', val);
    end
    
    %% Get/Set private methods
    function axes_size = get.axes_size(obj)
        axes_size = getpixelposition(obj.AxesHandle);
        axes_size = max(axes_size(3:4),2); % A size of 1, even if true, breaks things
    end
    
    function data_size = get.data_size(obj)
        if ~isempty(obj.Metadata)
            data_size = double([obj.Metadata{obj.current_frame}.ImageData.NumCols obj.Metadata{obj.current_frame}.ImageData.NumRows]);
        else % Before any file has been opened
            data_size = [1 1];
        end
    end
    
    function image_size = get.image_size(obj)
        % Displayed data is always 1 pixel less that figure size.  This is
        % a hack so that built-in pan works
        image_size=floor(obj.axes_size)-1;
    end
    
    function number_of_frames = get.number_of_frames(obj)
        number_of_frames = numel(obj.readerobj);
    end
    
    function total_toolbar_height = get.total_toolbar_height(obj)
        total_toolbar_height = obj.TOOLBAR_HEIGHT * (1 + (obj.number_of_frames>1));
    end
    
    %% Other public methods
    % Roughly the reverse of complexdatacursorfun.  Takes a lat/long and
    % places a datacursor-- as opposed to complexdatacursorfun which takes
    % a datacursor and computes a lat/long.
    function success = geojump(obj, lla)
        % Convert lat/long to row/column index
        row_col_pos = point_ground_to_slant(lla, obj.Metadata{obj.current_frame}).';
        % Check if valid within entire dataset
        if any(row_col_pos<1)||(row_col_pos(1)>obj.Metadata{obj.current_frame}.ImageData.NumRows)||...
                (row_col_pos(2)>obj.Metadata{obj.current_frame}.ImageData.NumCols)
            success = false;
            return; % Point is outside image
        else
            success = true;
        end
        % Check if point is within viewed area.  Otherwise move view.
        screenxlim=get(obj.AxesHandle,'Xlim');
        screenylim=get(obj.AxesHandle,'Ylim');
        coordUL=obj.axescoords2native([screenxlim(1) screenylim(1)]);
        coordLR=obj.axescoords2native([screenxlim(2) screenylim(2)]);
        if row_col_pos(1)<coordUL(2)||row_col_pos(1)>coordLR(2)||...
                row_col_pos(2)<coordUL(1)||row_col_pos(2)>coordLR(1)
            obj.CenterPos = fliplr(row_col_pos);
        end
        % Display new datatip
        newDatatip = createDatatip(datacursormode(obj.figure_handle), obj.image_handle);
        if verLessThan('matlab','8.4') % Unsupported feature that changed in MATLAB 2014b
            update(newDatatip, obj.native2axescoords(fliplr(row_col_pos)));
        else
            % In MATLAB 2014b, newDatatip.Cursor.Position = something
            % seemed to work.  In later versions, this functionality broke,
            % and the cursor would not always go where you set it, so the
            % following ugly hack was required to compute an offset between
            % where you asked it to go and where it really goes.
            desired_pos = round([obj.native2axescoords(fliplr(row_col_pos)) 0]);
            offset = [0 0 0];
            iters = 0;
            while any(abs(newDatatip.Cursor.Position - desired_pos) > 0.5) && iters<3
                newDatatip.Cursor.Position = desired_pos + offset;
                offset = offset + desired_pos - newDatatip.Cursor.Position;
                iters = iters + 1;
            end
        end
    end
    
    % Unlike setting public properties like Zoom or CenterPos, which set a
    % single parameter at a time, this function allows multiple viewing
    % parameters to be set simultaneously.
    function setView(obj, varargin)
        % Parse input parameters
        p = inputParser;
        p.KeepUnmatched=true;
        p.addParamValue('CenterPos', obj.CenterPos);
        p.addParamValue('Decimation', obj.Decimation);
        p.addParamValue('Frame', obj.Frame);
        p.addParamValue('Remap', obj.Remap);
        p.addParamValue('Zoom', obj.Zoom);
        p.parse(varargin{:});

        % Check input parameters
        if p.Results.Frame>obj.number_of_frames
            return;
        end
        obj.PreChangeViewFcn();
        if ~any(strcmp('Remap',p.UsingDefaults)) % Apply new remap only if requested
            obj.remaplist{p.Results.Frame}=p.Results.Remap;
        end
        if strncmpi(p.Results.Zoom,'fit',3)
            newzoom = max(ceil(obj.data_size./obj.axes_size));
        else
            newzoom = p.Results.Zoom;
        end
        % Load data
        if isequal(p.Results.CenterPos,obj.CenterPos)&&...
                isequal(p.Results.Decimation,obj.Decimation)&&...
                (newzoom==obj.Zoom)&&... % Frame is the only change
                ~isempty(obj.complex_image) % Not the first load of this file
            obj.loadNewDataNewFrameOnly(p.Results.Frame);
        else
            obj.current_frame=p.Results.Frame;
            obj.decimation_function = p.Results.Decimation;
            obj.loadNewData(newzoom, p.Results.CenterPos);
        end
        % Update comboBoxes
        if ~isempty(obj.decimateComboBox) % Without image processing toolbox, this won't exist
            if strcmpi(p.Results.Decimation,'none')
                obj.setSelectedItemNoCallback(obj.decimateComboBox, 'Decimate');
            else
                obj.setSelectedItemNoCallback(obj.decimateComboBox, p.Results.Decimation);
            end
        end
        obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
        if isa(p.Results.Remap, 'function_handle')
            remap_str = func2str(p.Results.Remap);
        else
            remap_str = p.Results.Remap;
        end
        obj.setSelectedItemNoCallback(obj.remapComboBox, remap_str);
        obj.setSelectedItemNoCallback(obj.zoomComboBox, num2str(newzoom));
    end
    
    % Convert coordinates in the MATLAB axes to coordinates in the full
    % image space
    function pos_out = axescoords2native(obj,pos_in)
        % Repeat a copy of the origin for each row in the input
        origin_matrix = repmat(obj.origin,size(pos_in,1),1);
        if(obj.scale<1)
            pos_out=pos_in+origin_matrix-1;
        else
            pos_out=(obj.scale*(pos_in-1))+origin_matrix;
        end
    end
    
    % Inverse of axescoords2native
    function pos_out = native2axescoords(obj,pos_in)
        % Repeat a copy of the origin for each row in the input
        origin_matrix = repmat(obj.origin,size(pos_in,1),1);
        if(obj.scale<1)
            pos_out=pos_in-origin_matrix+1;
        else
            pos_out=((pos_in-origin_matrix)/obj.scale)+1;
        end
    end
end

methods (Access = private)
    %% Layout manager for handling resize
    % Note: This resize function doesn't seem to work correctly if you
    % resize with the top bar only.
    function mitmResizeFcn(obj)
        % Manage layout for toolbars and axes within parent uipanel
        total_size = getpixelposition(obj.Parent);
        setpixelposition(obj.AxesHandle,...
            [1 1 total_size(3) max(total_size(4)-obj.total_toolbar_height,1)]);
        setpixelposition(obj.main_toolbar_hg,...
            [1 max(total_size(4)-obj.TOOLBAR_HEIGHT+1,1) total_size(3)+1 obj.TOOLBAR_HEIGHT]);
        obj.movie_toolbar.setVisible(obj.number_of_frames>1);
        setpixelposition(obj.movie_toolbar_hg,...
            [1 max(total_size(4)-obj.total_toolbar_height+1,1) total_size(3)+1 obj.TOOLBAR_HEIGHT]);
        obj.populateZoomComboBox(); % Keep all the powers of 2 up to widest zoom-out
        set(obj.AxesHandle,'XLim',[1 obj.axes_size(1)]); % Don't let figure automatically resize image
        set(obj.AxesHandle,'YLim',[1 obj.axes_size(2)]); % Keep at same scale
    end
    
    %% Navigation callbacks
    % Handle zoom in and zoom out
    function handleZoomClick(obj,src)
        if ~obj.panzoom
            obj.panzoom=true;
            zoomobj=zoom(src);
            if strcmp(get(zoomobj,'Direction'),'in')
                % Calculate requested scale.  Usually half the previous
                % zoom, but sometimes different with a zoom window
                screenxlim=get(obj.AxesHandle,'Xlim');
                screenylim=get(obj.AxesHandle,'Ylim');
                if (obj.scale>1) % Zooming in by loading new data
                    % Calculate requested scale.  Usually half the previous
                    % zoom, but sometimes different with a zoom window
                    corner1=obj.axescoords2native([screenxlim(1) screenylim(1)]);
                    corner2=obj.axescoords2native([screenxlim(2) screenylim(2)]);
                    newscale=max(round(max((corner2-corner1+1)./obj.axes_size)),1);
                else % Zooming in with MATLAB's built-in zoom, no new resolution or data
                    newscale=max(([diff(screenxlim) diff(screenylim)]+1)./obj.axes_size_prezoom);
                    obj.scale=2^round(log2(newscale));
                    obj.setSelectedItemNoCallback(obj.zoomComboBox, num2str(obj.scale));
                    obj.panzoom=false;
                    obj.PostChangeViewFcn();
                    return;
                end
            else
                newscale=obj.scale*2;
                if(abs(newscale-1)<eps)
                    newscale=1;
                end
                if(newscale<=1)
                    obj.scale=newscale;
                    obj.setSelectedItemNoCallback(obj.zoomComboBox, num2str(newscale));
                    obj.panzoom=false;
                    obj.PostChangeViewFcn();
                    return;
                end
            end
            obj.setSelectedItemNoCallback(obj.zoomComboBox, num2str(newscale));
            obj.loadNewData(newscale);
        end
    end
    
    % Built-in zoom-out doesn't generally work since only visible data is
    % available.  Only let built-in zoom-out work if we are zoomed in past
    % 1:1.  Otherwise handle zoom-out ourselves.
    function res = noZoomOut(obj)
        obj.PreChangeViewFcn(); % Run callback that should occur before each view change
        obj.axes_size_prezoom = obj.axes_size; % Save for handleZoomClick, will be used if obj.scale<=1
        zoomobj=zoom(obj.figure_handle);
        res=strcmp(get(zoomobj,'Direction'),'out')&&(obj.scale>=1);
        if res
            max_scale=max(ceil(obj.data_size./obj.axes_size));
            if(obj.scale==max_scale)
                return;
            end
            newscale=obj.scale*2;
            if(newscale>max_scale)
                newscale=max_scale;
            end
            obj.loadNewData(newscale);
            obj.setSelectedItemNoCallback(obj.zoomComboBox, num2str(obj.scale));
        end
    end
    
    % Handle panning
    function handlePanClick(obj)
        % Hand doesn't release current image until this whole funcion
        % completes.  So use a timer to run image update, so that hand
        % releases image before file read is performed.
        if ~obj.panzoom
            if(obj.scale>=1)
                obj.panzoom=true;
                start(obj.pan_timer)
            else
                obj.PostChangeViewFcn();
            end
        end
    end
    function handlePanClickTimer(obj)
        obj.loadNewData(obj.scale);
    end
    
    % We need one level of indirection here.  By assigning this
    % prePanAction function to the pan ActionPreCallback, rather than
    % PreChangeViewFcn directly, this allows the PreChangeViewFcn variable
    % to be evaluated at runtime, rather than assignment time, so a calling
    % function can change and edit this function.
    function prePanAction(obj)
        obj.PreChangeViewFcn();
    end
    
    % Makes the data cursor function display the coordinates in the
    % original data, not the subsampled data.  Also displays the value at
    % each pixel.
    function output_txt=complexdatacursorfun(obj,eventdata)
        screenpos=get(eventdata,'Position');
        pos=obj.axescoords2native(screenpos);
        value=double(squeeze(obj.complex_image{obj.current_frame}(screenpos(2),screenpos(1),:)).');
        output_txt={['X: ' num2str(round(pos(1))-1)],... % Zero-based to be consistent with SICD
            ['Y: ' num2str(round(pos(2))-1)]};
        pos_lla = point_slant_to_ground([pos(2); pos(1)], obj.Metadata{obj.current_frame});
        projectionDesc = 'projected to plane';
        if ~isempty(pos_lla)
            output_txt=[output_txt, {sprintf('Location (%s):',projectionDesc)}];
            output_txt=[output_txt, {sprintf('   Lat: %0.5f', pos_lla(1))}];
            output_txt=[output_txt, {sprintf('   Lon: %0.5f', pos_lla(2))}];
            output_txt=[output_txt, {sprintf('   Elev: %0.1f (HAE, %s)', pos_lla(3),'m')}];
        end
       
        output_txt=[output_txt, {['Value: ' num2str(value)]}];
        meta = obj.Metadata{obj.current_frame}; % For ease of notation
        if isfield(meta,'Grid')&&isfield(meta.Grid,'TimeCOAPoly')&&...
                ~isscalar(meta.Grid.TimeCOAPoly)
            coa_time = sicd_polyval2d(meta.Grid.TimeCOAPoly,pos(1),pos(2),meta);
            output_txt=[output_txt,{['COA Time: ' num2str(coa_time) ' (sec)']}];
        end
        % if isfield(meta,'RMA')&&isfield(meta.RMA,'INCA')&&...
        %         isfield(meta.RMA.INCA,'TimeCAPoly')
        %     coa_time = sicd_polyval2d(meta.RMA.INCA.TimeCAPoly(:).',pos(1),pos(2),meta);
        %     output_txt=[output_txt,{['CA Time: ' num2str(coa_time) ' (sec)']}];
        % end
    end

    % Read only the section of the file necessary for requested AOI and
    % display it.
    function loadNewData(obj, zoomscale, center)
        if nargin<3 % Use current center coords if none given
            center=obj.CenterPos;
        end
        obj.scale = zoomscale;

        currentxlim=[max(1,round(center(1)-(obj.scale*(obj.image_size(1)-1)/2))) 0];
        currentxlim(2)=currentxlim(1)+(obj.scale*obj.image_size(1))-1;
        if(currentxlim(2)>obj.data_size(1))
            currentxlim(2)=obj.data_size(1);
            currentxlim(1)=max(1,currentxlim(2)+1-(obj.scale*obj.image_size(1)));
        end
        currentylim=[max(1,round(center(2)-(obj.scale*(obj.image_size(2)-1)/2))) 0];
        currentylim(2)=currentylim(1)+(obj.scale*obj.image_size(2))-1;
        if(currentylim(2)>obj.data_size(2))
            currentylim(2)=obj.data_size(2);
            currentylim(1)=max(1,currentylim(2)+1-(obj.scale*obj.image_size(2)));
        end
        
        obj.origin=[currentxlim(1) currentylim(1)];
        obj.frames_loaded=false(size(obj.frames_loaded));
        obj.complex_image{obj.current_frame}=hg_mitm_viewer.reorientForDisplay(...
            obj.readerobj{obj.current_frame}.read_chip(currentxlim,currentylim,obj.scale*[1 1],obj.Decimation));
        obj.frames_loaded(obj.current_frame)=true;
        set(obj.image_handle,'CData',obj.makeDisplayable(...
            obj.complex_image{obj.current_frame},obj.remaplist{obj.current_frame}));
        for i=1:obj.number_of_frames
            obj.climlist{i}=[];
        end
        obj.setContrast();
        set(obj.AxesHandle,'XLim',[1 obj.axes_size(1)]);
        set(obj.AxesHandle,'YLim',[1 obj.axes_size(2)]);
        obj.PostChangeViewFcn();
        obj.panzoom=false;
        % Note that we cannot update any comboBoxes in this function.  If
        % we update the comboBoxes directly, we will trigger a callback
        % which will recursively call this function.  If we use the
        % setSelectedItemNoCallback function, then the temporary disabling
        % of the callback within that function will not work, if this
        % action was originally initiated from that callback.
    end
    
    % If we are only changing frame and not any other view parameters
    % (except remap), we cache the frames so that they can be pulled
    % quickly from memory, rather than read from file, so that movies can
    % be played.
    function loadNewDataNewFrameOnly(obj,frame_number_to_show)
        if frame_number_to_show>obj.number_of_frames
            return;
        end
        % Save previous frame info (colormap and color limits).  These may
        % be set externally outside of hg_mitm_viewer code, so we must
        % record them upon leaving a frame, rather than at set time.
        obj.cmaplist{obj.current_frame}=colormap(obj.AxesHandle);
        obj.climlist{obj.current_frame}=get(obj.AxesHandle,'CLim');
        % Show new figure
        obj.current_frame=frame_number_to_show;
        if ~obj.frames_loaded(obj.current_frame) % We must load from file
            currentxlim = obj.origin(1) + [0 min((obj.scale*obj.image_size(1)),obj.data_size(1))-1];
            currentylim = obj.origin(2) + [0 min((obj.scale*obj.image_size(2)),obj.data_size(2))-1];
            obj.complex_image{obj.current_frame}=hg_mitm_viewer.reorientForDisplay(...
                obj.readerobj{obj.current_frame}.read_chip(currentxlim,currentylim,obj.scale*[1 1],obj.Decimation));
            obj.frames_loaded(obj.current_frame)=true;
            set(obj.image_handle,'CData',obj.makeDisplayable(...
                obj.complex_image{obj.current_frame},obj.remaplist{obj.current_frame}));
            if isempty(obj.climlist{obj.current_frame})
                obj.setContrast();
            else
                obj.setContrast(obj.climlist{obj.current_frame});
            end
        else % If this frame was viewed before, just use data already in memory
            obj.applyNewRemapOnly(obj.remaplist{obj.current_frame});
        end
        % Its an unfortunate MATLAB limitation that each figure can only
        % have a single colormap.  Hopefully the parent figure isn't
        % already using another colormap in another axes.
        colormap(obj.AxesHandle, obj.cmaplist{obj.current_frame});
        obj.PostChangeViewFcn();
    end
    
    % If remap is the only thing we are changing in the display, we can
    % just use existing complex data in memory.  No need to fetch from
    % file.
    function applyNewRemapOnly(obj,new_remap)
        obj.remaplist{obj.current_frame}=new_remap;
        if ischar(obj.remaplist{obj.current_frame})
            obj.remaplist{obj.current_frame}=str2func(obj.remaplist{obj.current_frame});
        end
        if numel(obj.complex_image)>=obj.current_frame
            set(obj.image_handle,'CData',obj.makeDisplayable(...
                obj.complex_image{obj.current_frame},obj.remaplist{obj.current_frame}));
            obj.setContrast();
        end
    end

    % Do all processing to convert raw complex data into a real format
    % MATLAB can display
    function out = makeDisplayable(obj, in, remap)
        out = in;
        
        if isinteger(out) % Many remap functions won't work on complex ints
            out=single(out);
        end

        if ~isempty(obj.DataTransformFcn)
            out=feval(obj.DataTransformFcn,out);
        % For polarimetric data, provide a default color decomposition
        elseif isfield(obj.Metadata{obj.current_frame},'ImageFormation') && ...
                isfield(obj.Metadata{obj.current_frame}.ImageFormation,'TxRcvPolarizationProc') && ...
                iscell(obj.Metadata{obj.current_frame}.ImageFormation.TxRcvPolarizationProc) && ...
                (numel(obj.Metadata{obj.current_frame}.ImageFormation.TxRcvPolarizationProc)>1) && ...
                isempty(obj.DataTransformFcn)
            out=feval(str2func(obj.DEFAULT_POL_DECOMP),out,obj.Metadata{obj.current_frame});
        end
        
        if nargin>2&&~isempty(remap)
            % Apply remap per band
            try
                for i=1:size(out,3)
                    temp(:,:,i)=feval(remap,out(:,:,i));
                end
            catch
                error('HG_MITM_VIEWER:InvalidRemap',['Error applying remap function: ' func2str(remap)]);
            end
            out=temp; % Use intermediate variable to allow for datatype change
            if size(out,3)>1 && isfloat(out)% True color images must be between 0 and 1
                out=out-min(out(isfinite(out)));
                out=out/max(out(:));
            end
        end
    end
    
    %% Handle contrast
    % Scale between 0.5% and 99.5% of ordered values
    function new_clim=calcNewClim(obj)
        data_to_stretch=get(obj.image_handle,'CData');
        if isfloat(data_to_stretch)
            data_to_stretch=mean(data_to_stretch,3); % Normalize in multiband
            data_to_stretch(~isfinite(data_to_stretch))=0; % Remove Inf's and NaNs
            maxCData=max(data_to_stretch(:));
            minCData=min(data_to_stretch(:));
            if(~isempty(ver('images'))&&... % Image Processing toolbox required for stretchlim and mat2gray
                    isfinite(maxCData)&&isfinite(minCData))
                % At some point, I should probably write a version of this that
                % is Image Processing toolbox free, but its quick and easy this
                % way.
                gray_limits=stretchlim(mat2gray(data_to_stretch), .005);
                new_clim = gray_limits * double(maxCData - minCData);
                new_clim = new_clim + double(minCData);
            else
                datamean=mean(data_to_stretch(:));
                datastd=std(data_to_stretch(:));
                new_clim=[0 datamean+5*datastd];
            end
        elseif islogical(data_to_stretch)
            new_clim=[0 1];
        elseif isinteger(data_to_stretch)
            new_clim=[intmin(class(data_to_stretch)) intmax(class(data_to_stretch))];
        end
    end

    % Generalize contrast setting to truecolor data, not just grayscale
    function setContrast(obj, new_clims)
        if nargin<2
            new_clims=obj.calcNewClim();
        end
        olddata=get(obj.image_handle,'CData');
        istruecolor=(size(olddata,3)==3);
        if istruecolor&&~isinteger(olddata)
            olddata=min(new_clims(2),max(new_clims(1),olddata));
            set(obj.image_handle,'CData',(olddata-new_clims(1))/diff(new_clims));
        end
        set(obj.AxesHandle,'CLim',new_clims);
    end
    
    %% Support functions for dropboxes
    % Determine what items should be in the zoomComboBox.  Should be "Fit"
    % and every power of 2 under the "Fit" zoom level, down to 1.
    function populateZoomComboBox(obj)
        old_callback = get(handle(obj.zoomComboBox,'callbackproperties'),...
            'actionPerformedCallback');
        set(handle(obj.zoomComboBox,'callbackproperties'),...
            'actionPerformedCallback',[]); % Don't want this activity firing callback
        x=nextpow2(max(ceil(obj.data_size./obj.axes_size)));
        obj.zoomComboBox.removeAllItems();
        obj.zoomComboBox.addItem('Fit');
        while(x>0),
            x=x-1;
            obj.zoomComboBox.addItem(num2str(2^x));
        end;
        obj.zoomComboBox.setSelectedItem(num2str(obj.scale));
        % Turn ComboBox callback back on
        set(handle(obj.zoomComboBox,'callbackproperties'),...
            'actionPerformedCallback',old_callback);
    end
    
    % Handle any changes in the zoom comboBox
    function updateZoomComboBox(obj,eventdata)
        if ~isempty(obj.readerobj)
            obj.PreChangeViewFcn();
            if(get(get(eventdata,'Source'),'SelectedIndex')==0) % Fit
                newscale=max(ceil(obj.data_size./obj.axes_size));
                obj.loadNewData(newscale);
                set(get(eventdata,'Source'),'SelectedItem',num2str(obj.scale));
            else
                newscale=str2double(get(get(eventdata,'Source'),'SelectedItem'));
                if(newscale~=obj.scale)
                    if(newscale>0&&~rem(newscale,1)) % Check for valid value
                        obj.loadNewData(newscale);
                    else % Otherwise reset to previous value
                        set(get(eventdata,'Source'),'SelectedItem',num2str(obj.scale));
                    end
                end
            end
        end
    end

    % Handle any changes in the decimation comboBox
    function updateDecimateComboBox(obj,eventdata)
        if(get(get(eventdata,'Source'),'SelectedIndex')~=0)
            new_decimation=get(get(eventdata,'Source'),'SelectedItem');
        else
            new_decimation='none';
        end
        % If decimation has changed, redo decimation
        if(~strcmp(new_decimation,obj.decimation_function))
            obj.decimation_function=new_decimation;
            obj.loadNewData(obj.scale);
        end
    end

    % Determine what items should be in the frameComboBox.  Should be 1 to
    % the number of frames.
    function populateRemapComboBox(obj)
        remapFunctions = getremaplist();
        if isempty(remapFunctions)
            error('mitm_viewer:no_remap_found','No remap functions available.  Make sure function GETREMAPLIST is available.');
        end
        for i=1:length(remapFunctions)
            obj.remapComboBox.addItem(char(remapFunctions(i)));
        end
    end
        
    % Handle any changes in the remap comboBox
    function updateRemapComboBox(obj,eventdata)
        obj.applyNewRemapOnly(get(get(eventdata,'Source'),'SelectedItem'));
    end

    % Determine what items should be in the frameComboBox.  Should be 1 to
    % the number of frames.
    function populateFrameComboBox(obj)
        old_callback = get(handle(obj.frameComboBox,'callbackproperties'),...
            'actionPerformedCallback');
        set(handle(obj.frameComboBox,'callbackproperties'),...
            'actionPerformedCallback',[]); % Don't want this activity firing callback
        obj.frameComboBox.removeAllItems();
        for i = 1:obj.number_of_frames
            obj.frameComboBox.addItem(num2str(i));
        end;
        % Turn ComboBox callback back on
        set(handle(obj.frameComboBox,'callbackproperties'),...
            'actionPerformedCallback',old_callback);
    end
    
    % Handle any changes in the frame comboBox
    function updateFrameComboBox(obj,eventdata)
        obj.loadNewDataNewFrameOnly(get(get(eventdata,'Source'),'SelectedIndex')+1);
        obj.movie_toolbar.repaint;
        obj.movie_toolbar.revalidate;
    end

    %% Main toolbar functions
    % uitoolbar only allows the toolbar to be put in the top of a MATLAB
    % figure.  We want to use a toolbar in a subcomponent of a GUI, so we need
    % to build our own toolbar that reproduces the MATLAB figure toolbar, at
    % least for the buttons we require: zoom, pan, and data cursor.
    function add_main_toolbar_buttons(obj)
        obj.main_toolbar_buttons.current_tool = [];
        obj.main_toolbar_buttons.zoominButton = ...
            make_toolbar_button('tool_zoom_in.png', 'Zoom In');
        obj.main_toolbar.add(obj.main_toolbar_buttons.zoominButton);
        obj.main_toolbar_buttons.zoomoutButton = ...
            make_toolbar_button('tool_zoom_out.png', 'Zoom Out');
        obj.main_toolbar.add(obj.main_toolbar_buttons.zoomoutButton);
        obj.main_toolbar_buttons.handButton = ...
            make_toolbar_button('tool_hand.png', 'Pan');
        obj.main_toolbar.add(obj.main_toolbar_buttons.handButton);
        obj.main_toolbar_buttons.dataCursorButton = ...
            make_toolbar_button('tool_data_cursor.png', 'Data Cursor');
        obj.main_toolbar.add(obj.main_toolbar_buttons.dataCursorButton);
        
        % Basic function for creating all of our main toolbar toggle buttons.
        function java_obj = make_toolbar_button(icon_filename, tooltip)
            java_obj=javaObjectEDT('javax.swing.JToggleButton');
            java_obj.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
                fullfile(toolboxdir('matlab'),['icons' filesep icon_filename])));
            java_obj.setToolTipText(tooltip);
            java_obj.setFocusable(false); % Make consistent with MATLAB toolbars
            set(handle(java_obj,'callbackproperties'),'ActionPerformedCallback',...
                @(hbutton, eventStruct, hfig) obj.main_toolbar_callback(hbutton, eventStruct));
        end
    end

    % Patterned after putdowntext.m, the callback for the MATLAB figure
    % toolbar buttons as of 2012a
    function main_toolbar_callback(obj, hbutton, eventdata)
        if ~isempty(obj.main_toolbar_buttons.current_tool) &&...
                obj.main_toolbar_buttons.current_tool ~= eventdata.getSource
            obj.main_toolbar_buttons.current_tool.setSelected(0);
        end
        obj.main_toolbar_buttons.current_tool = eventdata.getSource;
        otherButtons = obj.main_toolbar.getComponents();
        for i =10:numel(otherButtons)
            try %if there are other defined toggle buttons in the toolbar outside of hg_mitm_viewer, shut em off.
                otherButtons(i).setSelected(0);
            end
        end
        if obj.main_toolbar_buttons.current_tool==obj.main_toolbar_buttons.zoominButton
            if hbutton.isSelected()
                zoom(obj.figure_handle,'inmode');
            else
                zoom(obj.figure_handle,'off')
            end
        elseif obj.main_toolbar_buttons.current_tool==obj.main_toolbar_buttons.zoomoutButton
            if hbutton.isSelected()
                zoom(obj.figure_handle,'outmode');
            else
                zoom(obj.figure_handle,'off')
            end
        elseif obj.main_toolbar_buttons.current_tool==obj.main_toolbar_buttons.handButton
            if hbutton.isSelected()
                pan(obj.figure_handle,'onkeepstyle')
            else
                pan(obj.figure_handle,'off');
            end
        elseif obj.main_toolbar_buttons.current_tool==obj.main_toolbar_buttons.dataCursorButton
            if hbutton.isSelected()
                datacursormode(obj.figure_handle,'on');
            else
                datacursormode(obj.figure_handle,'off');
            end
        end
    end

    %% Movie player toolbar functions
    % Create the movie toolbar
    % Derived from MPLAY on the MATLAB exchange
    function add_movie_toolbar_buttons(obj)
        % Load icons
        icons=load('mitm_icons');
        [icons.save,ignore,transp] = imread(fullfile(toolboxdir('matlab'), ...
            ['icons' filesep 'file_save.png']));
        icons.save = double(icons.save)/double(intmax(class(icons.save)));
        % Convert transparency to NaN style, consistent with mitm_icons.mat
        transp(transp~=0)=1; transp=double(transp); transp(transp==0)=NaN;
        for i=1:3, icons.save(:,:,i)=icons.save(:,:,i).*transp; end;

        % Draw buttons
        buttons.save = make_toolbar_button(icons.save, 'Save movie', ...
            @(hbutton, eventStruct, hfig) obj.movie_cb_save());
        obj.movie_toolbar.add(buttons.save);
        buttons.start = make_toolbar_button(icons.goto_start_default, 'Go to start', ...
            @(hbutton, eventStruct, hfig) obj.movie_cb_goto_start());
        obj.movie_toolbar.add(buttons.start);
        buttons.back = make_toolbar_button(icons.step_back, 'Step back', ...
            @(hbutton, eventStruct, hfig) obj.movie_cb_step_back());
        obj.movie_toolbar.add(buttons.back);
        buttons.play = make_toolbar_button(icons.play_on, 'Play/Pause', ...
            @(hbutton, eventStruct, hfig) obj.movie_cb_play());
        obj.movie_toolbar.add(buttons.play);
        buttons.forward = make_toolbar_button(icons.step_fwd, 'Step forward', ...
            @(hbutton, eventStruct, hfig) obj.movie_cb_step_fwd());
        obj.movie_toolbar.add(buttons.forward);
        buttons.end = make_toolbar_button(icons.goto_end_default, 'Go to end', ...
            @(hbutton, eventStruct, hfig) obj.movie_cb_goto_end());
        obj.movie_toolbar.add(buttons.end);
        buttons.repeat=javaObjectEDT('javax.swing.JToggleButton');
        buttons.repeat.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
            hg_mitm_viewer.im2javaWithAlpha(icons.loop_on)));
        buttons.repeat.setToolTipText('Repeat: On');
        buttons.repeat.setFocusable(false); % Make consistent with MATLAB toolbars
        set(handle(buttons.repeat,'callbackproperties'),'ActionPerformedCallback',...
            @(hbutton, eventStruct, hfig) obj.movie_cb_loop());
        obj.movie_toolbar.add(buttons.repeat);

        % Save inside toolbar appdata for later
        setappdata(obj.movie_toolbar_hg,'icons',icons);
        setappdata(obj.movie_toolbar_hg,'buttons',buttons);
        
        function java_obj = make_toolbar_button(icon, tooltip, callback)
            java_obj=javaObjectEDT('javax.swing.JButton');
            java_obj.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
                hg_mitm_viewer.im2javaWithAlpha(icon)));
            java_obj.setToolTipText(tooltip);
            java_obj.setFocusable(false); % Make consistent with MATLAB toolbars
            set(handle(java_obj,'callbackproperties'),'ActionPerformedCallback',callback);
        end
    end
    
    function timerTickFcn(obj)
        tb_app_data = getappdata(obj.movie_toolbar_hg);
        tb_app_data.buttons.play.setToolTipText('Pause');
        tb_app_data.buttons.play.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
                hg_mitm_viewer.im2javaWithAlpha(tb_app_data.icons.pause_default)));

        if(obj.current_frame<obj.number_of_frames)
            obj.loadNewDataNewFrameOnly(obj.current_frame+1);
            obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
        else
            obj.loadNewDataNewFrameOnly(1);
            obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
            % Is repeat playback turned off?
            if(tb_app_data.buttons.repeat.isSelected())
                stop(obj.movie_timer);
            end
        end
    end

    function timerStopFcn(obj)
        % Keep this here, not in movie_cb_stop
        % Could have stopped from stop button (eg, gone thru movie_cb_stop)
        % but also could have stopped here due to end of movie
        tb_app_data = getappdata(obj.movie_toolbar_hg);
        tb_app_data.buttons.play.setToolTipText('Resume');
        tb_app_data.buttons.play.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
            hg_mitm_viewer.im2javaWithAlpha(tb_app_data.icons.play_on)));
    end

        % Save movie button callback
    function movie_cb_save(obj)
        FPS = 4; % Currently a constant.  Perhaps in the future give users access via GUI
        isRunning = strcmp(get(obj.movie_timer,'Running'),'on');
        if ~isRunning
            filterspec={'*.avi','Audio Video Interleave (*.avi)';...
                        '*.gif','Animated GIF (*.gif)'};
            [filename,pathname,filterindex]=uiputfile(filterspec,'Save Movie As');
            if(filename)
                movietype = filterspec{filterindex,1}(3:end);
                if strcmp(movietype,'avi')
                    % VideoWriter was introduced in MATLAB R2010b.  This
                    % should be the way to write video files.  However, we
                    % want this code to be backward compatible to 2008b, so
                    % we use the old way.
                    aviobj=avifile([pathname filename], 'compression', 'None');
                    aviobj.fps=FPS;
                end
                % Make temporary figure just for saving movie
                imsize = size(obj.complex_image{obj.current_frame});
                temp_fig = figure('Position',[10 10 imsize([2 1])],'MenuBar','none','Toolbar','none');
                temp_ax = axes('Parent',temp_fig,'Position',[0 0 1 1]);
                temp_im = image(zeros(imsize),'Parent',temp_ax,'CDataMapping','scaled');
                set(temp_ax,'Visible','off','DataAspectRatio',[1 1 1]); % image changes these
                for frame_index=1:obj.number_of_frames
                    obj.loadNewDataNewFrameOnly(frame_index);
                    set(temp_im,'CData',get(obj.image_handle,'CData'));
                    set(temp_ax,'CLim',get(obj.AxesHandle,'CLim'));
                    set(temp_fig,'Colormap',get(obj.figure_handle,'Colormap'));
                    obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
                    try
                        fig_frame = getframe(temp_fig);
                        if strcmp(movietype,'avi')
                            aviobj = addframe(aviobj,fig_frame);
                        elseif strcmp(movietype,'gif')
                            if frame_index==1
                                imwrite(rgb2gray(frame2im(fig_frame)),...
                                    [pathname filename],'gif',...
                                    'DelayTime',1/FPS,'LoopCount',Inf);
                            else
                                imwrite(rgb2gray(frame2im(fig_frame)),...
                                    [pathname filename],'gif',...
                                    'WriteMode','append',...
                                    'DelayTime',1/FPS);
                            end
                        end
                    catch
                        le=lasterror;
                        if(strcmp(le.identifier,'MATLAB:capturescreen:RectangleMustBeAtLeastPartiallyOnScreen'))
                            warning('MITM_VIEWER:IMAGE_NOT_ON_CURRENT_SCREEN','MATLAB requires figure to be on the current screen.');
                            break;
                        else
                            rethrow(le);
                        end
                    end
                end
                if strcmp(movietype,'avi')
                    aviobj=close(aviobj);
                end
                close(temp_fig);
            end
        end
    end

    % Go to start button callback
    function movie_cb_goto_start(obj)
        if obj.current_frame~=1
            obj.loadNewDataNewFrameOnly(1);
            obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
        end
    end

    % Step one frame backward callback
    function movie_cb_step_back(obj)
        if obj.current_frame==1
            obj.loadNewDataNewFrameOnly(obj.number_of_frames);
        else
            obj.loadNewDataNewFrameOnly(obj.current_frame-1);
        end
        obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
    end

    % Play button callback
    function movie_cb_play(obj)
        isRunning = strcmp(get(obj.movie_timer,'Running'),'on');
        if isRunning
            stop(obj.movie_timer);
        else
            start(obj.movie_timer);
        end
    end

    % Step one frame forward callback
    function movie_cb_step_fwd(obj)
        if obj.current_frame==obj.number_of_frames;
            obj.loadNewDataNewFrameOnly(1);
        else
            obj.loadNewDataNewFrameOnly(obj.current_frame+1);
        end
        obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
    end

    function movie_cb_goto_end(obj)
    % Goto end button callback
        if obj.current_frame~=obj.number_of_frames
            obj.loadNewDataNewFrameOnly(obj.number_of_frames);
            obj.setSelectedIndexNoCallback(obj.frameComboBox, obj.current_frame-1);
        end
    end

    function movie_cb_loop(obj)
    % Repeat playback button callback
        tb_app_data = getappdata(obj.movie_toolbar_hg);
        if(tb_app_data.buttons.repeat.isSelected())
            tb_app_data.buttons.repeat.setToolTipText('Repeat: Off');
            tb_app_data.buttons.repeat.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
                hg_mitm_viewer.im2javaWithAlpha(tb_app_data.icons.loop_off)));

        else
            tb_app_data.buttons.repeat.setToolTipText('Repeat: On');
            tb_app_data.buttons.repeat.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
                hg_mitm_viewer.im2javaWithAlpha(tb_app_data.icons.loop_on)));
        end
    end
end

%% Helper functions
methods (Static, Access = private)
    % Reorients the data from the read_chip standard (1st dimension
    % aziumuth) so that MATLAB displays as energy-from-top
    function out = reorientForDisplay(in)
        permute_order=[2 1 3:ndims(in)];
        out = permute(in,permute_order);
    end
    
    % Derived from MATLAB's im2java, but modified to handle alpha channel
    function jimage = im2javaWithAlpha(rgb)
        mrows = size(rgb,1);
        ncols = size(rgb,2);
        alpha = 255*~isnan(rgb(:,:,1));
        packed = bitshift(uint32(alpha),24);
        rgb = rgb*255; % Input is required to be 8-bit, so scale
        packed = bitor(packed,bitshift(uint32(rgb(:,:,1)),16));
        packed = bitor(packed,bitshift(uint32(rgb(:,:,2)),8));
        packed = bitor(packed,uint32(rgb(:,:,3)));
        pixels = packed';
        mis = java.awt.image.MemoryImageSource(ncols,mrows,pixels(:),0,ncols);
        jimage = java.awt.Toolkit.getDefaultToolkit.createImage(mis);
    end
    
    % The following two functions set ComboBox selection without triggering
    % callback.  Generally all programmatic setting of the ComboBoxes
    % should be done through these functions.
    function setSelectedItemNoCallback(combobox, val)
        old_callback = get(handle(combobox,'callbackproperties'),...
            'actionPerformedCallback');
        set(handle(combobox,'callbackproperties'),...
            'actionPerformedCallback',[]);
        combobox.setSelectedItem(val);
        set(handle(combobox,'callbackproperties'),...
            'actionPerformedCallback',old_callback);
    end
    
    function setSelectedIndexNoCallback(combobox, val)
        old_callback = get(handle(combobox,'callbackproperties'),...
            'actionPerformedCallback');
        set(handle(combobox,'callbackproperties'),...
            'actionPerformedCallback',[]);
        combobox.setSelectedIndex(val);
        set(handle(combobox,'callbackproperties'),...
            'actionPerformedCallback',old_callback);
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
