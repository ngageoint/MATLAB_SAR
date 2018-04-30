function varargout = mitm_viewer( filename, varargin )
%MITM_VIEWER More-Image-Than-Memory complex image viewer
%   WARNING: This version of mitm_viewer is no longer maintained.  It is
%   only for version of MATLAB prior to 7.7 (2008b).  The new version of
%   mitm_viewer leverages the object oriented features of MATLAB introduced
%   in 7.6 (2008a) and javaObjectEDT which was introduced in 7.7 (2008b).
%
%   MITM_VIEWER (FILENAME, 'PropertyName', PropertyValue, ...) takes a
%   file of complex data (or set of files of identical sizes) and displays
%   them in a MATLAB IMSHOW/IMAGESC style viewer, even if the file is too
%   large to fit into memory.  Any file type handled by OPEN_READER can be
%   used.
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
%   The figure initially opens with an overview of the whole data file, but
%   the figure supports interactive zooming and panning.  Only as much
%   resolution and extent as are needed for the current figure size and
%   zoom level are read from file.  This allows for viewing of extremely
%   large datasets even on computers with limited memory.  Zooming can be
%   done with the magnifying glass tools in the toolbar or the scale factor
%   dropbox (also in the toolbar).  Single click zoom and window zoom
%   should both work.  Panning is done with the hand tool just as in a
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
%   decimated.)  Remap function take a single image as in input and output
%   a single image as an output.
%
%       Property name     Description
%       figureSize        Size of the figure to be drawn. If FIG_SIZE is
%                         not specified, the largest figure size that will
%                         fit in the current screen with an integer zoom
%                         factor is used.
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
%                         mode is 'navigate' or selectCallback is set, then
%                         this property is ignored.
%       undulationFilename   The name/path of the file containing WGS84 to
%                            geoid undulations.  This is used in the data
%                            cursor (see geoid_undulation) to display
%                            the height above ellipsoid rather then height
%                            above mean sea level.  If the file doesn't
%                            exist no HOE measurement is shown.
%
%   Dependencies: This function should run without errors with MATLAB only
%   and no toolboxes.  However, it has a few more features (optional
%   decimation and automatic colormap limiting) if you have the image
%   processing toolbox.  Also, the AOI selector will not work without the
%   image processing toolbox (althought it only requires
%   $matlabroot/toolbox/images/images/getrect.m and
%   $matalbroot/toolbox/images/images/private/getcurpt.m, so if you can
%   find them or equivalents, good for you.  I can't distribute them with
%   this code though because that would be illegal.)
%
%   Known issues:
%       Often locks up if you close figure while playing multiple frames.
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
p.FunctionName = 'MITM_VIEWER';
p.addParamValue('undulationFilename', 'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE');
p.parse(varargin{:});

%% Open files
persistent PATHNAME
if ((nargin<1)||isempty(filename)) % If no filename was give, use dialog box to ask for one
    try
        fileparts(PATHNAME); % Is PATHNAME valid
    catch
        PATHNAME=pwd;
    end
    [filename,pathname]=uigetfile(fullfile(PATHNAME,'*.*'),'MultiSelect','on');
    if(iscell(filename)), filename=sort(filename); end;
    PATHNAME=pathname;
else
    pathname=[];
end
if(iscell(filename)) % Multiple files requested
    number_of_filenames=length(filename);
    for i=1:length(filename)
        fullfilename{i}=[pathname filename{i}];
    end
elseif(filename)
    number_of_filenames=1;
    fullfilename{1}=[pathname filename];
else % filename=0.  Cancel was pressed, instead of a file being chosen.
    if ~strcmpi(p.Results.mode,'navigate')
        varargout={[],1};
    end
    return;
end
% Some files have multiple images within them.  Enumerate readers for all
% images in all files
readerobj={};
for i=1:number_of_filenames % Open readers for all files, and build cell array of all readers
    newreader=open_reader(fullfilename{i});
    readerobj=[readerobj newreader];
end
clear newreader;
% Some images are spread across multiple files
[reader_indices,readerobj]=group_reader_objs_by_pol(readerobj); % Consolidate reader objects into polarimetric sets
% Extract metadata for each image
metadata=cell(size(readerobj));
for i=1:length(readerobj)
    metadata{i}=readerobj{i}.get_meta();
    if((metadata{i}.ImageData.NumRows~=metadata{1}.ImageData.NumRows)||...
            (metadata{i}.ImageData.NumCols~=metadata{1}.ImageData.NumCols))
        warning('mitm_viewer:nonmatching_sizes','Sizes of all images passed to mitm_viewer are not the same.  Some features of viewer might not work.');
    end
end
datasize=double([metadata{1}.ImageData.NumCols metadata{1}.ImageData.NumRows]);
number_of_frames=length(readerobj);

%% Setup zoom scale and zoom scale dropbox
if isempty(p.Results.figureSize)
    screensize=get(0,'ScreenSize');
    fig_size=screensize(3:4)-screensize(1:2)-[20 100];  % Full screen size minus border
else
    fig_size=p.Results.figureSize;
end
scale=max(ceil(datasize./fig_size)); % Initial scale is as zoomed-out as will fit in figure
zoomComboBox=awtcreate('javax.swing.JComboBox');
awtinvoke(zoomComboBox,'setToolTipText','Zoom Level');
awtinvoke(zoomComboBox,'setMinimumSize',java.awt.Dimension(60,16));
awtinvoke(zoomComboBox,'setMaximumSize',java.awt.Dimension(60,1000));
awtinvoke(zoomComboBox,'setEditable',true);
repopulateZoomComboBox();
awtinvoke(zoomComboBox,'setSelectedItem',num2str(scale));
set(handle(zoomComboBox,'callbackproperties'),'actionPerformedCallback',@updateZoomComboBox);

%% Setup decimation types and dropbox
decimation_fun='none';
if(~isempty(ver('images'))) % Image Processing toolbox require for decimation feature in readers
    decimateComboBox=awtcreate('javax.swing.JComboBox');
    awtinvoke(decimateComboBox,'addItem','Decimate');
    awtinvoke(decimateComboBox,'addItem','max');
    awtinvoke(decimateComboBox,'addItem','mean');
    awtinvoke(decimateComboBox,'setToolTipText','Decimation Type');
    awtinvoke(decimateComboBox,'setMinimumSize',java.awt.Dimension(100,16));
    awtinvoke(decimateComboBox,'setMaximumSize',java.awt.Dimension(100,1000));
    awtinvoke(decimateComboBox,'setEditable',true);
    set(handle(decimateComboBox,'callbackproperties'),'ActionPerformedCallback',@updateDecimateComboBox);
end

%% Setup remap type and dropbox
remaplist=cell(number_of_frames,1); remaplist{1}=0;
remapFunctions = getremaplist();
remapComboBox=awtcreate('javax.swing.JComboBox');
for i=1:length(remapFunctions)
    name = char(remapFunctions(i));
    awtinvoke(remapComboBox,'addItem',char(remapFunctions(i)));
end
awtinvoke(remapComboBox,'setToolTipText','Remap Type');
awtinvoke(remapComboBox,'setMinimumSize',java.awt.Dimension(125,16));
awtinvoke(remapComboBox,'setMaximumSize',java.awt.Dimension(125,1000));
awtinvoke(remapComboBox,'setEditable',true);
if isempty(p.Results.initialRemap)
    if isempty(remapFunctions)
        error('mitm_viewer:no_remap_found','No remap functions available.  Make sure function GETREMAPLIST is available.');
    end    
    current_remap='densityremap'; % Default remap
else
    current_remap=p.Results.initialRemap;
end
if ischar(current_remap), current_remap=str2func(current_remap); end;
remap_index = find(strcmp(func2str(current_remap),remapFunctions));
if ~isempty(remap_index)
    awtinvoke(remapComboBox,'setSelectedIndex',remap_index-1);
else
    awtinvoke(remapComboBox,'setSelectedItem',func2str(current_remap));
end
set(handle(remapComboBox,'callbackproperties'),'ActionPerformedCallback',@updateRemapComboBox);

%% Setup frame selector
if(number_of_frames>1)
    frameComboBox=awtcreate('javax.swing.JComboBox');
    for i=1:number_of_frames
        awtinvoke(frameComboBox,'addItem',num2str(i));
    end
    awtinvoke(frameComboBox,'setSelectedIndex',p.Results.initialFrame-1);
    awtinvoke(frameComboBox,'setToolTipText','Frame Number');
    awtinvoke(frameComboBox,'setMinimumSize',java.awt.Dimension(50,16));
    awtinvoke(frameComboBox,'setMaximumSize',java.awt.Dimension(50,1000));
    set(handle(frameComboBox,'callbackproperties'),'ActionPerformedCallback',@updateFrameComboBox);
end

%% Metadata button
metaButton=awtcreate('javax.swing.JButton','Ljava.lang.String;', 'Metadata');
set(handle(metaButton,'callbackproperties'),'ActionPerformedCallback',@show_metadata);

%% Geojump button
geojumpButton=awtcreate('javax.swing.JButton','Ljava.lang.String;', 'Geojump');
set(handle(geojumpButton,'callbackproperties'),'ActionPerformedCallback',@geojump);

%% Return area (or points) of interest
if ~strcmpi(p.Results.mode,'navigate')
    switch lower(p.Results.mode)
        case 'aoi'
            button_str = 'Select AOI';
        case 'points'
            button_str = 'Select point(s)';
        case 'polygon'
            button_str = 'Select polygon';
    end
    aoiButton=awtcreate('javax.swing.JButton','Ljava.lang.String;', button_str);
    set(handle(aoiButton,'callbackproperties'),'ActionPerformedCallback',@select_AOI);
end

%% Setup initial figure
current_frame=p.Results.initialFrame;
frames_loaded=false(number_of_frames,1);
currentxlim=[];
currentylim=[];
complex_image{current_frame}=reorient_for_display(readerobj{current_frame}.read_chip(currentxlim,currentylim,[scale scale],decimation_fun));
frames_loaded(current_frame)=true;
origin=[1 1];
if isfield(metadata{current_frame},'CollectionInfo')&&isfield(metadata{current_frame}.CollectionInfo,'CoreName')
    figurename = metadata{current_frame}.CollectionInfo.CoreName;
else figurename = ['Image #' num2str(current_frame)];
end
fig_hand=figure('Name',figurename); % Figure handle
jftb=get(get(findall(fig_hand,'Type','uitoolbar'),'JavaContainer'),'ComponentPeer'); % Get Java handle for figure toolbar
if(number_of_frames>1)
    movietoolbar=create_movie_toolbar();
    movietimer = timer('ExecutionMode','fixedRate','TimerFcn', @TimerTickFcn, ...
        'StopFcn', @TimerStopFcn, 'BusyMode', 'drop', 'TasksToExecute', inf, 'Period', 0.5);
end
cmaplist=cell(number_of_frames,1);
climlist=cell(number_of_frames,1);
% im_hand=imshow(make_displayable(complex_image{current_frame},current_remap),[],'Border','tight'); % Image handle
% ax_hand=gca(); % Axes handle
% The next set of lines emulate imshow behavior.  We prefer to use imagesc
% to avoid the dependency on the Image Processing Toolbox.
%-----------------------
try
    im_hand=imagesc(make_displayable(complex_image{current_frame},current_remap));  % Image handle
catch
    error('mitm_viewer:invalid_remap',['Error applying remap function: ' func2str(current_remap)]);
end
cmaplist{1}=colormap(gray);
daspect([1 1 1]);
ax_hand=gca(); % Axes handle
climlist{1}=calc_new_clim;
set_contrast(climlist{1});
set(ax_hand,'Position',[0 0 1 1]);
set(ax_hand,'TickDir','out');
%----------------------
for i=1:number_of_frames, cmaplist{i}=cmaplist{1}; end;
set(fig_hand,'ResizeFcn',@myresizefcn); % Requires axes handle before it can be added
fig_size=size(complex_image{current_frame})+1; % Make figure 1 pixel larger than image so that built-in MATLAB pan works
fig_size=fig_size([2 1]);
set(fig_hand,'Position',[20 20 fig_size]);
set(fig_hand,'DeleteFcn',@cleanmemory); % Java/MATLAB memory leak if you don't explicitly clear components
% WARNING: Adding Java components to the MATLAB toolbar is an unsupported
% feature.  This may not work in a future version.
if(~isempty(jftb))
    jftb.addSeparator;
    jftb.add(zoomComboBox);
    if(~isempty(ver('images'))), jftb.add(decimateComboBox); end;
    jftb.add(remapComboBox);
    jftb.addSeparator;
    jftb.add(metaButton);
    try
        if ~isempty(point_slant_to_ground([1; 1], metadata{current_frame}, 'projectToDEM', false))
            jftb.add(geojumpButton);
        end
    end
    if(~strcmpi(p.Results.mode,'navigate'))
        jftb.addSeparator;
        jftb.add(aoiButton);
    end
    jftb.repaint;
    jftb.revalidate;
end
if(number_of_frames>1)
    jmtb=[];
    while(isempty(jmtb)) % Takes a while to add new toolbar.  Why?  I don't know...
        jmtb=get(get(findall(fig_hand,'tag','movietoolbar'),'JavaContainer'),'ComponentPeer'); % Get Java handle for movie toolbar
        drawnow();
    end
    if(~isempty(jmtb))
        jmtb.addSeparator;
        jmtb.add(frameComboBox);
        jmtb.repaint;
        jmtb.revalidate;
    end
end

%% Setup figure toolbar callbacks
set(datacursormode(fig_hand),'UpdateFcn',@complexdatacursorfun);
set(zoom(fig_hand),'ActionPostCallback',@handleZoomClick);
set(zoom(fig_hand),'ButtonDownFilter',@noZoomOut);
set(pan(fig_hand),'ActionPostCallback',@handlePanClick);
zooming=false; panning=false;
panTimer=timer('TimerFcn',@handlePanClickTimer,'StartDelay',.1);
% Will selection parameters be returned in output argument?  If so, wait
if(~strcmpi(p.Results.mode,'navigate')&&isempty(p.Results.selectCallback))
    % Should disable metadata button here until AOI selection is finished
    waitfor(fig_hand,'UserData','AOI_finished');
    if(~exist('varargout','var'))
        varargout={[],1};
    end
end


%% Callbacks for figure toolbar button functions
    % Makes the data cursor function display the coordinates in the
    % original data, not the subsampled data.  Also displays the value at
    % each pixel.
    function output_txt=complexdatacursorfun(obj,eventdata)
        screenpos=get(eventdata,'Position');
        pos=screen2native(screenpos);
        value=double(squeeze(complex_image{current_frame}(screenpos(2),screenpos(1),:)).');
        output_txt={['X: ' num2str(round(pos(1)))],...
            ['Y: ' num2str(round(pos(2)))]};
        meta_struct=metadata{current_frame};
        
%         Currently disabled because it makes data cursor much slower.
%         try
%             % Try to project the point to a DEM (most accurate projection)
%             pos_lla = point_slant_to_ground([pos(2); pos(1)], meta_struct, 'projectToDEM', true);
%             projectionDesc = 'projected to DEM';
%         catch
            % Failing the projection to a DEM, just fall back to projecting
            % to a plane.
            pos_lla = point_slant_to_ground([pos(2); pos(1)], meta_struct, 'projectToDEM', false);
            projectionDesc = 'projected to plane';
%         end
        if ~isempty(pos_lla)
            output_txt={output_txt{:}, sprintf('Location (%s):',projectionDesc)};
            output_txt={output_txt{:}, sprintf('   Lat: %0.5f', pos_lla(1))};
            output_txt={output_txt{:}, sprintf('   Lon: %0.5f', pos_lla(2))};
            output_txt={output_txt{:}, sprintf('   Elev: %0.1f (HAE, %s)', pos_lla(3),'m')};
            if (exist(p.Results.undulationFilename, 'file'))
              undulation = geoid_undulation(pos_lla(1), pos_lla(2), 'undulationFilename', p.Results.undulationFilename);
              elev_msl = pos_lla(3) - undulation;
              output_txt={output_txt{:}, sprintf('   Elev: %.1f (MSL, %s)', elev_msl,'m')};
            end
        end
       
        if isreal(complex_image{current_frame})
            output_txt={output_txt{:},['Value: ' num2str(value)]};
        else
            output_txt={output_txt{:},...
            ['Complex: ' num2str(value)],...
            ['Abs: ' num2str(abs(value))],...
            ['Angle: ' num2str(angle(value))]};
        end
    end

    % Roughly the reverse of complexdatacursorfun.  Takes a lat/long and
    % places a datacursor-- as opposed to complexdatacursorfun which takes
    % a datacursor and computes a lat/long.
    function geojump(obj,eventdata)
        meta = metadata{current_frame}; % For easier referencing
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
        % Convert lat/long to row/column index
        row_col_pos = point_ground_to_slant(lla, meta, 'projectToDEM', false).';
        % Check if valid within entire dataset
        if any(row_col_pos<1)||(row_col_pos(1)>meta.ImageData.NumRows)||...
                (row_col_pos(2)>meta.ImageData.NumCols)
            msgbox('Specified point is not in image!');
            return; % Point is outside image
        end
        % Check if point is within viewed area.  Otherwise move view.
        screenxlim=get(ax_hand,'Xlim');
        screenylim=get(ax_hand,'Ylim');
        coordUL=screen2native([screenxlim(1) screenylim(1)]);
        coordLR=screen2native([screenxlim(2) screenylim(2)]);
        if row_col_pos(1)<coordUL(2)||row_col_pos(1)>coordLR(2)||...
            row_col_pos(2)<coordUL(1)||row_col_pos(2)>coordLR(1)
            readImageAndRedisplay(scale, fliplr(row_col_pos));
        end
        % Display new datatip
        update(createDatatip(datacursormode(fig_hand), im_hand),...
            native2screen(fliplr(row_col_pos)));
        
        % Opposite of screen2native, but only local since we only use it here
        function pos_out = native2screen(pos_in)
            if(scale<1)
                pos_out=pos_in-origin+1;
            else
                pos_out=((pos_in-origin)/scale)+1;
            end
        end
    end

    % Handle zoom in and zoom out
    function handleZoomClick(obj,eventdata)
        if ~(zooming||panning)
            zooming=true;
            zoomobj=zoom(obj);
            if strcmp(get(zoomobj,'Direction'),'in')
                % Calculate requested scale.  Usually half the previous
                % zoom, but sometimes different with a zoom window
                screenxlim=get(ax_hand,'Xlim');
                screenylim=get(ax_hand,'Ylim');
                if (scale>1) % Zooming in by loading new data
                    % Calculate requested scale.  Usually half the previous
                    % zoom, but sometimes different with a zoom window
                    corner1=screen2native([screenxlim(1) screenylim(1)]);
                    corner2=screen2native([screenxlim(2) screenylim(2)]);
                    newscale=max(round(max((corner2-corner1+1)./fig_size)),1);
                    % We call readImageAndRedisplay indirectly through ComboBox
                    % update callback so ComboBox will be updated before file
                    % is read.  Since MATLAB built-in zoom will show low-res
                    % zoom almost immediately, and file reads take a little
                    % time, we want ComboBox to be up-to-date with what is
                    % showing on figure.
                    awtinvoke(zoomComboBox,'setSelectedItem',num2str(newscale));
                    % Since the previous set function should trigger
                    % readImageAndRedisplay through the ComboBox update
                    % callback, so no need to do it again at the end of this
                    % function; just return.
                    return;
                else % Zooming in with MATLAB's built-in zoom, no new resolution or data
                    newscale=max([diff(screenylim) diff(screenxlim)]./size(complex_image{current_frame}));
                    scale=2^round(log2(newscale));
                    awtinvoke(zoomComboBox,'setSelectedItem',num2str(newscale));
                    zooming=false;
                    return;
                end
            else
                newscale=scale*2;
                if(abs(newscale-1)<eps)
                    newscale=1;
                end
                if(newscale<=1)
                    scale=newscale;
                    awtinvoke(zoomComboBox,'setSelectedItem',num2str(newscale));
                    zooming=false;
                    return;
                end
            end
            readImageAndRedisplay(newscale);
        end
    end

    % Handle panning
    function handlePanClick(obj,eventdata)
        % Hand doesn't release current image until this whole funcion
        % completes.  So use a timer to run image update, so that hand
        % releases image before file read is performed.
        if ~(zooming||panning)
            panning=true;
            if(scale>=1), start(panTimer); end
        end
    end
    function handlePanClickTimer(obj,eventdata)
        readImageAndRedisplay(scale);
    end

    % For a given zoom level, read only the section of the image necessary
    % and display it.
    function readImageAndRedisplay(zoomscale, center)
        if nargin<2
            center=get_center_coords();
        end
        scale=zoomscale;
        % Make displayed data 1 pixel less that figure size.  This is a
        % hack so that built-in pan works
        usable_fig_size=fig_size-1;
        
        currentxlim=[max(1,center(1)-(scale*floor(usable_fig_size(1)/2))) 0];
        currentxlim(2)=currentxlim(1)+(scale*usable_fig_size(1))-1;
        if(currentxlim(2)>datasize(1))
            currentxlim(2)=datasize(1);
            currentxlim(1)=max(1,currentxlim(2)+1-(scale*usable_fig_size(1)));
        end
        currentylim=[max(1,center(2)-(scale*floor(usable_fig_size(2)/2))) 0];
        currentylim(2)=currentylim(1)+(scale*usable_fig_size(2))-1;
        if(currentylim(2)>datasize(2))
            currentylim(2)=datasize(2);
            currentylim(1)=max(1,currentylim(2)+1-(scale*usable_fig_size(2)));
        end
        
        origin=[currentxlim(1) currentylim(1)];
        frames_loaded=false(size(frames_loaded));
        complex_image{current_frame}=reorient_for_display(readerobj{current_frame}.read_chip(currentxlim,currentylim,[scale scale],decimation_fun));
        frames_loaded(current_frame)=true;
        set(im_hand,'CData',make_displayable(complex_image{current_frame},current_remap));
        for i=1:number_of_frames
            climlist{i}=[];
        end
        set_contrast();
        set(ax_hand,'XLim',[1 fig_size(1)]);
        set(ax_hand,'YLim',[1 fig_size(2)]);
        awtinvoke(zoomComboBox,'setSelectedItem',num2str(scale));
        zooming=false; panning=false;
    end

    % Built-in zoom-out doesn't generally work since only visible data is
    % available.  Only let built-in zoom-out work if we are zoomed in past
    % 1:1.  Otherwise handle zoom-out ourselves.
    function res = noZoomOut(obj,eventdata)
        zoomobj=zoom(fig_hand);
        max_scale=max(ceil(datasize./fig_size));
        res=strcmp(get(zoomobj,'Direction'),'out')&&(scale>=1);
        if (res)
            if(scale==max_scale)
                return;
            end
            newscale=scale*2;
            if(newscale>max_scale)
                newscale=max_scale;
            end
            readImageAndRedisplay(newscale);
        end
    end

%% Callbacks for dropboxes 
    % Handle any changes in the zoom comboBox
    function updateZoomComboBox(obj,eventdata)
        if(get(get(eventdata,'Source'),'SelectedIndex')==0)
            % Fit overview into window
            readImageAndRedisplay(max(ceil(datasize./fig_size)));
        else
            newscale=str2double(get(get(eventdata,'Source'),'SelectedItem'));
            if(newscale~=scale)
                if(newscale>0&&~rem(newscale,1))
                    readImageAndRedisplay(newscale);
                else
                    set(get(eventdata,'Source'),'SelectedItem',num2str(scale));
                end
            end
        end
    end

    % Handle any changes in the decimation comboBox
    function updateDecimateComboBox(obj,eventdata)
        if(get(get(eventdata,'Source'),'SelectedIndex')~=0)
            new_decimation_fun=get(get(eventdata,'Source'),'SelectedItem');
        else
            new_decimation_fun='none';
        end
        % If decimation has changed, redo decimation
        if(~strcmp(new_decimation_fun,decimation_fun))
            decimation_fun=new_decimation_fun;
            readImageAndRedisplay(scale);
            set_contrast();
        end
    end

    % Handle any changes in the remap comboBox
    function updateRemapComboBox(obj,eventdata)
        current_remap=get(get(eventdata,'Source'),'SelectedItem');
        if ischar(current_remap), current_remap=str2func(current_remap); end;
        if exist('im_hand','var') % Initial set sometimes happens before image is created
            set(im_hand,'CData',make_displayable(complex_image{current_frame},current_remap));
            set_contrast();
        end
    end

    % Handle any changes in the frame comboBox
    function updateFrameComboBox(obj,eventdata)
        if(get(eventdata,'modifiers')) % Changed from the GUI, and not programatically
            showNewFrame(get(get(eventdata,'Source'),'SelectedIndex')+1);
        end
        jmtb.repaint;
        jmtb.revalidate;
    end

    % Display image metadata
    function show_metadata(obj,eventdata)
        % Just opens the MATLAB variable editor for the metadata.  Assumes
        % user is at the MATLAB commandline with no debugger.
        % assignin('base','metadata',metadata{current_frame});
        % openvar('metadata');
        % Brings up new figure with tree view of metadata
        % Doesn't handle cell arrays or display 2-D matrices as well as
        % the MATLAB variable editor, but cleaner because its standalone
        % and doesn't write to base MATLAB workspace.
        MetaViewer(metadata{current_frame});
    end

    function select_AOI(obj,eventdatda)
        % Get AOI
        if strcmpi(p.Results.mode,'aoi') % Rectangular AOI
            aoi_screen=[0 0 0 0];
            while(any(aoi_screen(3:4)==0)),
                aoi_screen=getrect(ax_hand);
            end
        elseif strcmpi(p.Results.mode,'polygon')
            [x,y]=getline(ax_hand);
        elseif (p.Results.numPoints>0)&&isfinite(p.Results.numPoints) % Predetermined number of points
            [x,y]=ginput(p.Results.numPoints);
        else % Indefinite number of points
            [x,y]=getpts(ax_hand);
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
            % Convert to image coordinate space
            if strcmpi(p.Results.mode,'aoi')
                ulcorner=round(screen2native(aoi_screen(1:2)));
                lrcorner=ulcorner+round(aoi_screen(3:4)*scale)-1;
                % Constrain to be within image dimensions
                ulcorner=max(1,ulcorner);
                lrcorner=min(datasize,lrcorner);
                % Return values (UL corner and size) to calling function
                varargout{1}(1:2)=ulcorner;
                varargout{1}(3:4)=lrcorner-ulcorner+1;
            else % Series of points
                for i=1:length(x)
                    varargout{1}(i, 1:2)=round(screen2native([x(i) y(i)]));
                end
            end
            varargout{2}=reader_indices{current_frame}; % Which frame was AOI selected in?
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
                    set(fig_hand,'UserData','AOI_finished'); % If you want to keep window, do this
                    drawnow; % If you don't do this, everything freezes when you disable button
                    awtinvoke(aoiButton,'setEnabled',false);
                end
            end
        else
            delete(h); if exist('h2','var'), delete(h2); end;
        end
    end

    function showNewFrame(frame_number_to_show)
        % Save previous frame info (colormap, color limits, and remap)
        cmaplist{current_frame}=colormap();
        climlist{current_frame}=get(ax_hand,'CLim');
        if(get(remapComboBox,'SelectedIndex')>0)
            remaplist{current_frame}=get(remapComboBox,'SelectedIndex');
        else
            remaplist{current_frame}=get(remapComboBox,'SelectedItem');
        end
        % Show new figure
        current_frame=frame_number_to_show;
        if(~frames_loaded(current_frame))
            complex_image{current_frame}=reorient_for_display(readerobj{current_frame}.read_chip(currentxlim,currentylim,[scale scale],decimation_fun));
            frames_loaded(current_frame)=true;
            set(im_hand,'CData',make_displayable(complex_image{current_frame},current_remap));
            if isempty(climlist{current_frame})
                set_contrast();
            else
                set_contrast(climlist{current_frame});
            end
        else % If this frame was viewed before, recall display info
            set(im_hand,'CData',make_displayable(complex_image{current_frame},current_remap));
            set_contrast(climlist{current_frame});
            if isnumeric(remaplist{current_frame})
                awtinvoke(remapComboBox,'setSelectedIndex',remaplist{current_frame});
            else
                awtinvoke(remapComboBox,'setSelectedItem',remaplist{current_frame});
            end
        end
        colormap(cmaplist{current_frame});
        if isfield(metadata{current_frame},'CollectionInfo')&&isfield(metadata{current_frame}.CollectionInfo,'CoreName')
            figurename = metadata{current_frame}.CollectionInfo.CoreName;
        else figurename = ['Image #' num2str(current_frame)];
        end
        set(fig_hand,'Name',figurename); % Figure handle
    end

%% Movie player toolbar functions
    % Create the movie toolbar
    % Taken from MPLAY on the MATLAB exchange
    function tb_hand=create_movie_toolbar()
        tb_hand=uitoolbar('Parent',fig_hand,'tag','movietoolbar'); % Draw second toolbar
        % Load icons
        icons=load('mitm_icons');
        setappdata(tb_hand,'icons',icons); % Store icons in toolbar appdata

        % Draw buttons
        uipushtool(tb_hand, 'cdata',...
            getfield(load(fullfile(matlabroot, 'toolbox\matlab\icons\savedoc.mat'),'cdata'),'cdata'),...
            'tooltip','Save movie', 'click', @cb_save_movie);
        uipushtool(tb_hand, 'cdata', icons.goto_start_default, ...
            'tooltip','Go to start', 'click', @cb_goto_start);
        uipushtool(tb_hand, 'cdata', icons.step_back, ...
            'tooltip','Step back', 'click', @cb_step_back);
        uipushtool(tb_hand, 'cdata', icons.play_on, ...
            'tooltip','Play', 'tag','Play/Pause', 'click', @cb_play);
        uipushtool(tb_hand, 'cdata', icons.step_fwd, ...
            'tooltip','Step forward', 'click', @cb_step_fwd);
        uipushtool(tb_hand, 'cdata', icons.goto_end_default, ...
            'tooltip','Go to end', 'click', @cb_goto_end);
        uitoggletool(tb_hand, 'cdata', icons.loop_on, ...
            'tooltip','Repeat: On', 'tag','loopbutton', ...
            'state', 'on', 'click', @cb_loop);
    end

    function TimerTickFcn(hco, user)
        hPlay = findobj(movietoolbar, 'tag','Play/Pause');
        tb_app_data = getappdata(movietoolbar);
        set(hPlay, 'tooltip', 'Resume', 'cdata', tb_app_data.icons.pause_default);

        if(current_frame<number_of_frames)
            showNewFrame(current_frame+1);
            awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
        else
            showNewFrame(1);
            awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
            % Is repeat playback turned off?
            if(~strcmp(get(findobj(movietoolbar, 'tag','loopbutton'),'State'),'on'))
                stop(movietimer);
            end
        end
    end

    function TimerStopFcn(hco, user)
        % Keep this here, not in cb_stop
        % Could have stopped from stop button (eg, gone thru cb_stop)
        % but also could have stopped here due to end of movie
        hPlay = findobj(movietoolbar, 'tag','Play/Pause');
        tb_app_data = getappdata(movietoolbar);
        set(hPlay, 'tooltip', 'Resume', 'cdata', tb_app_data.icons.play_on);
    end

    % Save movie button callback
    function cb_save_movie(hbutton, eventStruct, hfig)
        FPS = 4; % Currently a constant.  Perhaps in the future give users access via GUI
        isRunning = strcmp(get(movietimer,'Running'),'on');
        if ~isRunning
            filterspec={'*.avi','Audio Video Interleave (*.avi)';...
                        '*.gif','Animated GIF (*.gif)'};
            [filename,pathname,filterindex]=uiputfile(filterspec,'Save Movie As');
            if(filename)
                movietype = filterspec{filterindex,1}(3:end);
                if strcmp(movietype,'avi')
                    aviobj=avifile([pathname filename], 'compression', 'None');
                    aviobj.fps=FPS;
                end
                for frame_index=1:number_of_frames
                    showNewFrame(frame_index);
                    awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
                    try
                        fig_frame=getframe(fig_hand);
                        if strcmp(movietype,'avi')
                            aviobj =addframe(aviobj,fig_frame);
                        elseif strcmp(movietype,'gif')
                            if frame_index==1
                                imwrite(rgb2gray(frame2im(fig_frame)),...
                                    [pathname filename],'gif',...
                                    'DelayTime',1/FPS,'LoopCount',Inf);
                            else
                                imwrite(rgb2gray(frame2im(fig_frame)),...
                                    [pathname filename],'gif',...
                                    'WriteMode','append');
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
            end
        end
    end

    % Go to start button callback
    function cb_goto_start(hbutton, eventStruct, hfig)
        if current_frame~=1
            showNewFrame(1);
            awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
        end
    end

    % Step one frame backward callback
    function cb_step_back(hbutton, eventStruct, hfig)
        if current_frame==1
            showNewFrame(number_of_frames);
        else
            showNewFrame(current_frame-1);
        end
        awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
    end

    % Play button callback
    function cb_play(hbutton, eventStruct, hfig)
        isRunning = strcmp(get(movietimer,'Running'),'on');
        if isRunning
            stop(movietimer);
        else
            start(movietimer);
        end
    end

    % Step one frame forward callback
    function cb_step_fwd(hbutton, eventStruct, hfig)
        if current_frame==number_of_frames;
            showNewFrame(1);
        else
            showNewFrame(current_frame+1);
        end
        awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
    end

    function cb_goto_end(hbutton, eventStruct, hfig)
    % Goto end button callback
        if current_frame~=number_of_frames
            showNewFrame(number_of_frames);
            awtinvoke(frameComboBox,'setSelectedIndex',current_frame-1);
        end
    end

    function cb_loop(hbutton, eventStruct, hfig)
    % Repeat playback button callback
        tb_app_data = getappdata(movietoolbar);
        if(strcmp(get(hbutton,'State'),'off'))
            set(hbutton, 'TooltipString', 'Repeat: Off',...
                'cdata', tb_app_data.icons.loop_off);
        else
            set(hbutton, 'TooltipString', 'Repeat: On',...
                'cdata', tb_app_data.icons.loop_on);
        end
    end

%% Callbacks for figure object
    % Handle figure resizing.
    function myresizefcn(obj,eventdata)
        fig_size=get(fig_hand,'Position'); fig_size=fig_size(3:4); % Make sure figure size is up to date
        repopulateZoomComboBox(); % Keep all the powers of 2 up to widest zoom-out
        set(ax_hand,'XLim',[1 fig_size(1)]); % Don't let figure automatically resize image
        set(ax_hand,'YLim',[1 fig_size(2)]); % Keep at same scale
    end

    % Figure close callback
    function cleanmemory(obj, eventdata)
        % Java/MATLAB memory leak can occure if you don't explicitly clear components
        cb1=[];
        cb1h=[];
        cb2=[];
        cb2h=[];
        % Clear timers
        while(strcmp(get(panTimer,'Running'),'on'));
            stop(panTimer);
        end
        delete(panTimer);
        if(number_of_frames>1)
            while(strcmp(get(movietimer,'Running'),'on'));
                stop(movietimer);
            end
            delete(movietimer);
        end
        % Close readers
        for i=1:number_of_frames
            readerobj{i}.close();
        end
        clear readerobj; % Releases memory mapped files
    end

%% Helper functions
    % Convert coordinates in the figure to coordinates in the full image
    % space
    function pos_out = screen2native(pos_in)
        if(scale<1)
            pos_out=pos_in+origin-1;
        else
            pos_out=(scale*(pos_in-1))+origin;
        end
    end

    % Get center of current view in coordinates from the full image
    function coords = get_center_coords()
        screenxlim=get(ax_hand,'Xlim');
        screenylim=get(ax_hand,'Ylim');
        coords=round(screen2native([mean(screenxlim) mean(screenylim)]));
    end

    % Determine what items should be in the zoom ComboBox.  Should be "Fit"
    % and every power of 2 under the "Fit" zoom level, down to 1.
    function repopulateZoomComboBox()
        x=nextpow2(max(ceil(datasize./fig_size)));
        awtinvoke(zoomComboBox,'removeAllItems');
        awtinvoke(zoomComboBox,'addItem','Fit');
        while(x>0),
            x=x-1;
            awtinvoke(zoomComboBox,'addItem',num2str(2^x));
        end;
        awtinvoke(zoomComboBox,'setSelectedItem',num2str(scale));
    end

    % Scale between 0.5% and 99.5% of ordered values
    function new_clim=calc_new_clim()
        data_to_stretch=get(im_hand,'CData');
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
    function set_contrast(new_clims)
        if nargin<1
            new_clims=calc_new_clim;
        end
        olddata=get(im_hand,'CData');
        istruecolor=(size(olddata,3)==3);
        if istruecolor&&~isinteger(olddata)
            olddata=min(new_clims(2),max(new_clims(1),olddata));
            set(im_hand,'CData',(olddata-new_clims(1))/diff(new_clims));
        end
        set(ax_hand,'CLim',new_clims);
    end

    % Reorients the data from the read_chip standard (1st dimension
    % aziumuth) so that imshow/imagesc displays as energy-from-top
    function out = reorient_for_display(in)
        permute_order=[2 1 3:ndims(in)];
        out = permute(in,permute_order);
    end

    % Do all processing to convert raw complex data into a format MATLAB
    % can display
    function out = make_displayable(in, remap, decomposition)
        out = in;
        
        if isinteger(out) % Many remap functions won't work on complex ints
            out=single(out);
        end
        
        if nargin<3
            switch size(out,3)
                case 1 % Single band image; nothing to do
                case 2 % Dual-pol (quasi-Pauli)
                   co_index=[find(strcmpi('H:H',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc))...
                      find(strcmpi('V:V',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc))];
                   if length(co_index)>1 % Catch HH/VV
                      cross_index=co_index(2);
                      co_index=co_index(1);
                   else
                      cross_index=[find(strcmpi('H:V',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc))...
                        find(strcmpi('V:H',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc))];
                   end
                   out=cat(3,abs(out(:,:,co_index)),...
                      abs(out(:,:,cross_index)),...
                      abs(out(:,:,co_index)));
                case 3 % RGB image; nothing to do
                case 4 % Quad-pol: default to Pauli decomposition
                    HH_index=find(strcmpi('H:H',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc));
                    HV_index=find(strcmpi('H:V',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc));
                    VH_index=find(strcmpi('V:H',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc));
                    VV_index=find(strcmpi('V:V',metadata{current_frame}.ImageFormation.TxRcvPolarizationProc));
                    out=cat(3,abs(out(:,:,HH_index)-out(:,:,VV_index)),...
                        abs(out(:,:,HV_index)+out(:,:,VH_index))/2,...
                        abs(out(:,:,HH_index)+out(:,:,VV_index)));
            end
        else
            out=feval(decomposition,out);
        end
        
        if nargin>1&&~isempty(remap)
            % Apply remap per band
            for i=1:size(out,3)
                temp(:,:,i)=feval(remap,out(:,:,i));
            end
            out=temp; % Use intermediate variable to allow for datatype change
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////