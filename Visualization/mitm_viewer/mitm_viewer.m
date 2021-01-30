function varargout = mitm_viewer( filename, varargin )
%MITM_VIEWER More-Image-Than-Memory complex image viewer
%   MITM_VIEWER (FILENAME, 'PropertyName', PropertyValue, ...) takes a
%   file of complex data (or set of files of identical sizes) and displays
%   them in a MATLAB IMSHOW/IMAGESC style viewer, even if the file is too
%   large to fit into memory.  Any file type handled by OPEN_READER can be
%   used.
%
%   Really just a convenient GUI wrapper for the HG_MITM_VIEWER component.
%
%       Property name     Description
%       figureSize        Size of the figure to be drawn (if a new figure
%                         is to be created).  If figureSize is not
%                         specified, the largest figure size that will fit
%                         in the current screen with an integer zoom factor
%                         is used.  This propert is ignored if the parent
%                         property is set.
%       initialFrame      Set the initial frame viewed.  Default is 1.
%       initialRemap      Set the initial remap function.  Default is
%                         'densityremap'.
%       mode              'navigate' (default), 'aoi', 'points', or
%                         'polygon'. Navigate just browses images.  'aoi'
%                         adds a button that allows the users to select a
%                         rectangular AOI which is returned as the output
%                         parameter.  'points' selects a single point or
%                         set of points in the image.  'polygon' actually
%                         behaves the same as 'points' with an indefinite
%                         number points, but draws lines between the
%                         points.
%       numPoints         If mode is 'points', this is the number of points
%                         to select.  Default is 1.  Setting to 0 or Inf
%                         will allow an indefinite number of points (same
%                         use as MATLAB getpts function).  If mode is not
%                         'points', this parameter is ignored.
%       selectCallback    Function handle for function that should be run
%                         after selection.  Callback is a function that
%                         takes a single argument, the selection
%                         parameters.  If selectCallback is set, multiple
%                         AOI/points can be chosen.  If not, only one
%                         AOI/point is chosen and returned in the output
%                         argument.
%       closeAfterSelect  Default is false.  Determines whether the figure
%                         should close after the AOI/point is chosen.  If
%                         mode is 'navigate', or if parent or
%                         selectCallback are set, then this property is
%                         ignored.
%       undulationFilename   The name/path of the file containing WGS84 to
%                            geoid undulations.  This is used in the data
%                            cursor (see geoid_undulation) to display
%                            the height above ellipsoid rather then height
%                            above mean sea level.  If the file doesn't
%                            exist no HOE measurement is shown.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('figureSize', [], @(x) isequal(size(x),[1 2])||isequal(size(x),[2 1]));
p.addParamValue('mode', 'navigate', @(x) any(strcmpi(x,{'navigate', 'aoi', 'points','polygon'})));
p.addParamValue('numPoints',1,@(x) isscalar(x));
p.addParamValue('initialFrame',1, @(x) isscalar(x)&&isnumeric(x)&&(x>0));
p.addParamValue('initialRemap',[], @(x) exist(x)>0);
p.addParamValue('closeAfterSelect', false);
p.addParamValue('selectCallback',[],@(x) isempty(x)||isa(x,'function_handle'))
p.addParamValue('undulationFilename', 'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE');
p.FunctionName = 'MITM_VIEWER';
p.parse(varargin{:});

%% Open files
if ((nargin<1)||isempty(filename)) % If no filename was give, use dialog box to ask for one
    % Recall last interactively selected path used
    if ispref('matlab_sar_toolbox','last_used_directory')
        pathname = getpref('matlab_sar_toolbox','last_used_directory');
        if ~ischar(pathname)||~exist(pathname,'dir')
            pathname = pwd;
        end
    else
        pathname = pwd;
    end
    [filename,pathname]=uigetfile(sar_file_extensions('complex'),...
        'Select Input File', pathname,'MultiSelect','on');
    if(iscell(filename)), filename=sort(filename); end;
    setpref('matlab_sar_toolbox','last_used_directory',pathname);
else
    pathname=[];
end
if(isinmem(filename)) % Filename is actually an array in memory
    fullfilename{1}=filename;
elseif(iscell(filename)) % Multiple files requested
    for i=1:length(filename)
        fullfilename{i}=[pathname filename{i}];
    end
elseif(filename) % Single filename
    fullfilename{1}=[pathname filename];
else % filename=0.  Cancel was pressed, instead of a file being chosen.
    if ~strcmpi(p.Results.mode,'navigate')
        varargout={[],1};
    end
    return;
end

%% Setup layout
fig_hand=figure('MenuBar','none','Toolbar','none');
uip_hand = uipanel('Parent',fig_hand,'Position',[0 0 1 1],'BorderType','none');
mitm_hand = hg_mitm_viewer(uip_hand);
% Save image button
saveButton=javaObjectEDT('javax.swing.JButton');
saveButton.setIcon(javaObjectEDT('javax.swing.ImageIcon',...
    fullfile(toolboxdir('matlab'),['icons' filesep 'file_save.png'])));
saveButton.setToolTipText('Save As...');
saveButton.setFocusable(false); % Make consistent with MATLAB toolbars
set(handle(saveButton,'callbackproperties'),'ActionPerformedCallback',...
    @(hbutton, eventStruct, hfig) save_image());
% Metadata button
metaButton=javaObjectEDT('javax.swing.JButton','Metadata');
set(handle(metaButton,'callbackproperties'),'ActionPerformedCallback',@(obj,eventdata) show_metadata());
% Geojump button
geojumpButton=javaObjectEDT('javax.swing.JButton','Geojump');
set(handle(geojumpButton,'callbackproperties'),'ActionPerformedCallback',@(obj,eventdata) geojump());
% AOI button
if ~strcmpi(p.Results.mode,'navigate')
    switch lower(p.Results.mode)
        case 'aoi'
            button_str = 'Select AOI';
        case 'points'
            button_str = 'Select point(s)';
        case 'polygon'
            button_str = 'Select polygon';
    end
    aoiButton=javaObjectEDT('javax.swing.JButton',button_str);
    set(handle(aoiButton,'callbackproperties'),'ActionPerformedCallback',@(obj,eventdata) select_AOI());
end

% Open file
mitm_hand.Frame = p.Results.initialFrame;
mitm_hand.openFile(fullfilename, true); 
metadata = mitm_hand.Metadata;

% Add buttons to toolbar
main_toolbar = mitm_hand.main_toolbar;
main_toolbar.addSeparator;
main_toolbar.add(metaButton);
try
    if ~isempty(point_slant_to_ground([1; 1], metadata{mitm_hand.Frame}))
        main_toolbar.add(geojumpButton);
    end
end
if(~strcmpi(p.Results.mode,'navigate'))
    main_toolbar.addSeparator;
    main_toolbar.add(aoiButton);
end
main_toolbar.addSeparator;
main_toolbar.add(saveButton);
main_toolbar.repaint;
main_toolbar.revalidate;

%% Set figure size appropriate for data
datasize=double([metadata{mitm_hand.Frame}.ImageData.NumCols metadata{mitm_hand.Frame}.ImageData.NumRows]);
if numel(metadata)==1 % Ugly.  Hardcoded to match the toolbar height found in hg_mitm_viewer.
    toolbar_height = 25;
else
    toolbar_height = 50;
end
if isempty(p.Results.figureSize)
    screensize=get(0,'ScreenSize');
    % Image must fit in full screen size minus toolbars and border
    max_img_size=((screensize(3:4)-screensize(1:2))-[0 toolbar_height]) - [20 50];
else
    max_img_size=p.Results.figureSize;
end
if isfield(metadata{mitm_hand.Frame},'CollectionInfo')&&...
        isfield(metadata{mitm_hand.Frame}.CollectionInfo,'CoreName')
    fig_name_str = metadata{mitm_hand.Frame}.CollectionInfo.CoreName;
elseif (numel(metadata)>1)
    fig_name_str = ['Image #' num2str(mitm_hand.Frame)];
else
    fig_name_str = '';
end
scale = max(ceil(datasize./max_img_size)); % Initial scale is lowest integer decimation that will fit
% Calculation of desired size is very specific to internal workings of hg_mitm_viewer
new_size = ceil(datasize./scale) + [0 toolbar_height] + 1;
set(fig_hand,'Name',fig_name_str,'Position',[20 20 new_size]);
% The resize callback executes in parallel with the current code.  Make
% sure resize completes, before we continue, since the following zoom is
% dependent on the new size.
tid = tic;
while(~isequal(getpixelposition(fig_hand),[20 20 new_size]))&&(toc(tid)<1)
end
mitm_hand.Zoom = 'fit';
if ~isempty(p.Results.initialRemap)
    mitm_hand.Remap = p.Results.initialRemap;
end

%% Will selection parameters be returned in output argument?  If so, wait
if(~strcmpi(p.Results.mode,'navigate')&&isempty(p.Results.selectCallback))
    % Should disable metadata button here until AOI selection is finished
    waitfor(uip_hand,'UserData','AOI_finished');
    if(~exist('varargout','var'))
        varargout={[],1};
    end
end


%% Callbacks for figure toolbar button functions
    % All of the work here is temporarily removing the toolbar from the
    % figure, so it doesn't get saved.  We assume that we must have a
    % figure with nothing else other than the mitm_viewer.
    function save_image()
        old_resize = get(uip_hand,'ResizeFcn');
        old_visible = get(mitm_hand.movie_toolbar,'Visible');
        set(uip_hand,'ResizeFcn',[]);
        fig_pos = getpixelposition(fig_hand);
        axes_pos = getpixelposition(mitm_hand.AxesHandle);
        setpixelposition(fig_hand,[fig_pos(1:2) axes_pos(3:4)]);            
        setpixelposition(mitm_hand.AxesHandle,axes_pos);            
        filemenufcn(fig_hand,'FileSaveAs');
        set(mitm_hand.movie_toolbar,'Visible',old_visible); % Why does FileSaveAs change this?
        set(uip_hand,'ResizeFcn',old_resize);
        setpixelposition(fig_hand,fig_pos); % Replace toolbar
    end

    % Roughly the reverse of complexdatacursorfun.  Takes a lat/long and
    % places a datacursor-- as opposed to complexdatacursorfun which takes
    % a datacursor and computes a lat/long.
    function geojump()
        full_meta = mitm_hand.Metadata;
        meta = full_meta{mitm_hand.Frame}; % For easier referencing
        % Get lat/long from user
        if isfield(meta,'GeoData')&&isfield(meta.GeoData,'SCP')
            % Default elevation is elevation of scene center
            def_el = num2str(meta.GeoData.SCP.LLH.HAE);
        else
            def_el = '';
        end
        answers = inputdlg({'Latitude','Longitude','Elevation (HAE, m)'},'Geojump',1,{'','',def_el});
        if isempty(answers), return, end; % Cancel was pressed
        lla = [latlonnum(answers{1});...
            latlonnum(answers{2});...
            str2double(answers{3})];
        if any(isnan(lla))
            msgbox('Invalid Lat/Lon format!');
            return;
        end
        if ~mitm_hand.geojump(lla)
            msgbox('Specified point is not in image!');
        end
    end

    % Display image metadata
    function show_metadata
        % Just opens the MATLAB variable editor for the metadata.  Assumes
        % user is at the MATLAB commandline with no debugger.
        % assignin('base','metadata',metadata{mitm_hand.Frame});
        % openvar('metadata');
        % Brings up new figure with tree view of metadata
        % Doesn't handle cell arrays or display 2-D matrices as well as
        % the MATLAB variable editor, but cleaner because its standalone
        % and doesn't write to base MATLAB workspace.
        MetaViewer(metadata{mitm_hand.Frame});
    end

    function select_AOI
        % Get AOI
        if strcmpi(p.Results.mode,'aoi') % Rectangular AOI
            aoi_screen=[0 0 0 0];
            while(any(aoi_screen(3:4)==0)),
                aoi_screen=getrect(mitm_hand.AxesHandle);
            end
        elseif strcmpi(p.Results.mode,'polygon')
            [x,y]=getline(mitm_hand.AxesHandle);
        elseif (p.Results.numPoints>0)&&isfinite(p.Results.numPoints) % Predetermined number of points
            [x,y]=ginput(p.Results.numPoints);
        else % Indefinite number of points
            [x,y]=getpts(mitm_hand.AxesHandle);
        end
        % Draw AOI
        if strcmpi(p.Results.mode,'aoi')
            h=rectangle('Position',aoi_screen,'LineWidth',2,'EdgeColor','r');
        elseif strcmpi(p.Results.mode,'polygon')
            h=line([x; x(1)],[y; y(1)],zeros(size(x)+[1 0]),'Color', 'r', 'LineWidth', 2);
        else % Points
            % Emulate getpts-style markers
            h=line(x,y,zeros(size(x)),'Color', 'c', 'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 9);
            h2=line(x,y,'Color', 'm', 'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 9);
        end
        % Confirm selection
        if ~strcmpi(p.Results.mode,'aoi')
            if p.Results.numPoints~=1
                question_string='Are these the points you wish to select?';
            else
                question_string='Is this the point you wish to select?';
            end
        else
            question_string='Is this the area of interest you wish to select?';
        end
        button=questdlg(question_string,'Confirm AOI','Yes','No','Yes');
        if(strcmp(button,'Yes'))
            varargout={[],1};
            % Convert to image coordinate space
            if strcmpi(p.Results.mode,'aoi')
                ulcorner=round(mitm_hand.axescoords2native(aoi_screen(1:2)));
                lrcorner=round(mitm_hand.axescoords2native(aoi_screen(1:2) + aoi_screen(3:4) - 1));
                % Constrain to be within image dimensions
                ulcorner=max(1,ulcorner);
                lrcorner=min(datasize,lrcorner);
                % Return values (UL corner and size) to calling function
                varargout{1}(1:2)=ulcorner;
                varargout{1}(3:4)=lrcorner-ulcorner+1;
            else % Series of points
                for i=1:length(x)
                    varargout{1}(i, 1:2)=round(mitm_hand.axescoords2native([x(i) y(i)]));
                end
            end
            varargout{2}=mitm_hand.FrameSubFileIndices{mitm_hand.Frame}; % Which frame was AOI selected in?
            if isa(p.Results.selectCallback,'function_handle')
                try
                    p.Results.selectCallback(varargout{:});
                catch
                    delete(h); if exist('h2','var'), delete(h2); end;
                    rethrow(lasterror);
                end
                delete(h); if exist('h2','var'), delete(h2); end;
            else
                if (p.Results.closeAfterSelect)
                    close(fig_hand);
                else                
                    set(uip_hand,'UserData','AOI_finished'); % If you want to keep window, do this
                    drawnow; % If you don't do this, everything freezes when you disable button
                    aoiButton.setEnabled(false);
                end
            end
        else
            delete(h); if exist('h2','var'), delete(h2); end;
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////