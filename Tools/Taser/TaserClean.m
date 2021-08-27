function TaserClean( filename )
%TASERCLEAN Taser SAR Data Analyis GUI (minimilistic version)
%
% Minimilistic version of TASER interface.  TASER itself is not a
% processing algorithm, but rather is a GUI that allows access to a number
% of different MATLAB SAR processing algorithms and visualization routines.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup layout
taser_str = 'Taser (Clean)';
fig_hand=figure('Name', taser_str, 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'Toolbar', 'none', 'DockControls', 'off');
initial_pos = getpixelposition(fig_hand);
setpixelposition(fig_hand,[initial_pos(1), initial_pos(2)+initial_pos(4), [250 1]]);
uip_hand = uipanel('Parent',fig_hand,'Position',[0 0 1 1],'BorderType','none');
mitm_hand = hg_mitm_viewer(uip_hand); % For complex data
set(mitm_hand.main_toolbar,'Visible',false);
phd = struct('uip_children', [], 'filename', '', 'meta', []); % For phase history

% Get tool list
algo_fid = fopen([fileparts(mfilename('fullpath')) filesep 'Algorithms.txt']);
algo_cell = textscan(algo_fid,'%s %s %s %s %s','delimiter',',','HeaderLines',1);
algorithm_info = cell2struct([algo_cell{:}],{'TextName','CallName','DataType','Selected','Polarimetric'},2);
fclose(algo_fid);

% Menubar
file_uimenu = uimenu(fig_hand, 'Label', 'File');
uimenu(file_uimenu, 'Label', 'Open...', ...
    'Callback', @(obj, event) open_file());
snapshot_menu = uimenu(file_uimenu, 'Label', 'Save snapshot...', 'Enable', 'off', ...
    'Callback', @(obj, event) save_image());
meta_menu = uimenu(file_uimenu, 'Label', 'Metadata', 'Separator', 'on');
metaicon_menu = uimenu(meta_menu, 'Label', 'Metaicon','Enable', 'off', ...
    'Callback',@(obj, event) metaicon_callback());
metatree_menu = uimenu(meta_menu, 'Label', 'Tree view','Enable', 'off', ...
   'Callback',@(obj, event) metaviewer_callback()); 
uimenu(file_uimenu, 'Label', 'Exit', 'Separator', 'on', ...
    'Callback',@(obj, event) close(fig_hand));
navigate_uimenu = uimenu(fig_hand, 'Label', 'Navigate');
geojump_menu = uimenu(navigate_uimenu, 'Label', 'Geojump...', 'Enable', 'off', ...
    'Callback', @(obj, event) geojump_gui());
seg_map_menu = uimenu(navigate_uimenu, 'Label', 'Segment Map...', ...
   'Callback', @(obj, event) plot_segment_map(),'Enable', 'off'); 
algorithms_uimenu = uimenu(fig_hand, 'Label', 'Tools');
populate_algorithm_menu();
help_uimenu = uimenu(fig_hand, 'Label', 'Help');
uimenu(help_uimenu, 'Label', 'About', 'Callback', @(obj, event) taser_about());

% Open file if passed in
if nargin>0
    open_file(filename);
end

    %% Callbacks
    function open_file(filename)
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
            [filename,pathname]=uigetfile(sar_file_extensions({'complex','phd'}),...
                'Open SAR Data File',pathname,'MultiSelect','on');
            if(iscell(filename)), filename=sort(filename); end;
            setpref('matlab_sar_toolbox','last_used_directory',pathname);
        else % Path was already passed in with filename
            pathname='';
        end
        fullfilename = {}; % Reset if files were previously open
        if(iscell(filename)) % Multiple files requested
            for j=1:length(filename)
                fullfilename{j}=[pathname filename{j}];
            end
        elseif(filename)
            fullfilename{1}=[pathname filename];
        else % filename=0.  Cancel was pressed, instead of a file being chosen.
            return;
        end
        
        % Close previous data
        mitm_hand.close(); % If old data was complex
        delete(phd.uip_children); % If old data was phase history
        phd = struct('uip_children', [], 'filename', '', 'meta', []);

        % Open new data
        if ~isempty(guess_ph_format(fullfilename{1})) % Phase history data
            old_pos = getpixelposition(fig_hand);
            NEW_SIZE = [600 400];
            setpixelposition(fig_hand,[old_pos(1), ...
                old_pos(2)+old_pos(4)-NEW_SIZE(2), NEW_SIZE]);
            set([metaicon_menu metatree_menu], 'Enable', 'on');
            set([snapshot_menu geojump_menu seg_map_menu], 'Enable', 'off');
            set(mitm_hand.main_toolbar,'Visible',false);
            phd.filename = fullfilename{1}; % Save for later
            % Open file
            ph_reader_object = open_ph_reader(phd.filename);
            populate_algorithm_menu();
            phd.meta = ph_reader_object.get_meta();
            NumPulses = phd.meta.Data.Channel(1).NumVectors;
            if isfield(phd.meta,'CollectionID')&&...
                    isfield(phd.meta.CollectionID,'CoreName')
                corename = phd.meta.CollectionID.CoreName;
            else
                corename = '';
            end
            fig_name_str = [taser_str ': ' corename];
            % TODO: Get these from preferences
            NumMPSDPulses = 100;
            ProcessingPulses = 100;
            
            % MPSD
            [freq mpsddata] = MPSD(ph_reader_object, false, 1, ...
                NumPulses, floor(NumPulses/NumMPSDPulses), ProcessingPulses);
            phd.uip_children(1) = axes('Parent',uip_hand,'Position',[.1 .1 .85 .35]);
            plot(phd.uip_children(1),freq./10^9,mpsddata);
            xlabel(phd.uip_children(1),'Frequency (GHz)');
            ylabel(phd.uip_children(1),'Power (dB)');
            xlim(phd.uip_children(1),[min(freq./10^9) max(freq./10^9)]);
            title(phd.uip_children(1),['Mean Power Spectral Density for : ' corename]);
            
            % Meta icon
            phd.uip_children(2) = axes('Parent',uip_hand,'Position',[0 .55 1 .4]);
            MetaIcon_PHD(phd.meta, 'handle', phd.uip_children(2));
            ph_reader_object.close();
        else % Complex data
            % set figure size appropriate for data
            mitm_hand.openFile(fullfilename, true);
            set(mitm_hand.main_toolbar,'Visible',true);
            metadata = mitm_hand.Metadata;
            populate_algorithm_menu();
            set([snapshot_menu metaicon_menu metatree_menu geojump_menu], 'Enable', 'on');
            if length(metadata)>1
                set(seg_map_menu,'Enable', 'on');
            else
                set(seg_map_menu,'Enable', 'off');
            end
            
            datasize=double([metadata{mitm_hand.Frame}.ImageData.NumCols metadata{mitm_hand.Frame}.ImageData.NumRows]);
            if numel(metadata)==1 % Ugly.  Hardcoded to match the toolbar height found in hg_mitm_viewer.
                toolbar_height = 25;
            else
                toolbar_height = 50;
            end
            screensize=get(0,'ScreenSize');
            % Image must fit in full screen size minus toolbars and border
            max_img_size=((screensize(3:4)-screensize(1:2))-[0 toolbar_height]) - [20 50];
            if isfield(metadata{mitm_hand.Frame},'CollectionInfo')&&...
                    isfield(metadata{mitm_hand.Frame}.CollectionInfo,'CoreName')
                fig_name_str = [taser_str ': ' metadata{mitm_hand.Frame}.CollectionInfo.CoreName];
            else
                fig_name_str = taser_str;
            end
            % Initial scale is lowest integer decimation that will fit
            scale = max(ceil(datasize./max_img_size));
            % Calculation of desired size is very specific to internal workings of hg_mitm_viewer
            new_size = ceil(datasize./scale) + [0 toolbar_height] + 1;
            setpixelposition(fig_hand,[20 20 new_size]);
            % The resize callback executes in parallel with the current code.
            % Make sure resize completes, before we continue, since the
            % following zoom is dependent on the new size.
            tid = tic;
            while(~isequal(getpixelposition(fig_hand),[20 20 new_size]))&&(toc(tid)<1)
            end
            mitm_hand.Zoom = 'fit';
        end
        set(fig_hand,'Name',fig_name_str);
    end

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

    function metaicon_callback()
        fig1 = figure('MenuBar', 'none', 'Toolbar', 'none');
        pos1 = getpixelposition(fig1);
        pos1(3:4) = 200*[1 1];
        setpixelposition(fig1,pos1);
        ax1 = axes('Parent',fig1,'Position',[0 0 1 1]);
        if isempty(phd.meta)
            MetaIcon_Complex(mitm_hand.Metadata{mitm_hand.Frame},'handle',ax1);
        else
            MetaIcon_PHD(phd.filename,'handle',ax1);
        end
    end

    function metaviewer_callback()
        if isempty(phd.meta)
            MetaViewer(mitm_hand.Metadata{mitm_hand.Frame});
        else
            MetaViewer(phd.meta);
        end
    end

    function geojump_gui()
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

    function plot_segment_map()
       figure;
       hold on;
       for j = 1:length(mitm_hand.Metadata)
          meta = mitm_hand.Metadata{j};
          ICP = meta.GeoData.ImageCorners.ICP;
          Lats = [ICP.FRFC.Lat ICP.FRLC.Lat ICP.LRLC.Lat ICP.LRFC.Lat ICP.FRFC.Lat];
          Lons = [ICP.FRFC.Lon ICP.FRLC.Lon ICP.LRLC.Lon ICP.LRFC.Lon ICP.FRFC.Lon];
          text(mean(Lons(1:4)),mean(Lats(1:4)),num2str(j));
          plot(Lons,Lats);
       end
       hold off;

       xlabel('Longitude (deg)');
       ylabel('Latitude (deg)');
       title(['Segment Map for Collect: ' meta.CollectionInfo.CoreName(1:16)]);
    end

    function run_tool(toolinfo)
        if ~isempty(phd.filename) % Phase history data loaded.  No AOI.
            feval(toolinfo.CallName, 'filename', phd.filename);
        elseif ~isempty(mitm_hand.Metadata) % Complex data loaded.  Visible data is used as AOI.
            % Check for multichannel data
            filenames = mitm_hand.Filename(mitm_hand.FrameFileIndices{mitm_hand.Frame});
            segments = mitm_hand.FrameSubFileIndices{mitm_hand.Frame};
            if numel(filenames)>1||numel(segments)>1 % Data is multi-pol
                if strcmpi(toolinfo.Polarimetric,'no')
                    [choice,ok] = listdlg('SelectionMode','single',...
                        'PromptString',{'Algorithm only processes', 'one data channel.', 'Process single channel?'}, ...
                        'ListString',mitm_hand.Metadata{mitm_hand.Frame}.ImageFormation.TxRcvPolarizationProc,...
                        'Name','Warning','ListSize',[160 80]);
                    if ~ok, return; end;
                    filenames = filenames{choice};
                    segments = segments(choice);
                else
                    filenames = unique(filenames);
                end
            else
                filenames = filenames{1}; % Simple string more likely to be handled than cell array
            end
            feval(toolinfo.CallName, 'filename', filenames, ... % Run tool
                'aoi', mitm_hand.Position, 'segment', segments);
        else % No data loaded
            eval(toolinfo.CallName);
        end
    end

    function populate_algorithm_menu()
        delete(get(algorithms_uimenu,'Children'));
        
        data_is_phd = ~isempty(phd.filename);
        data_is_complex = ~isempty(mitm_hand.Metadata);
        data_is_multipol = data_is_complex && ...
            (numel(mitm_hand.FrameFileIndices{mitm_hand.Frame})>1 ||...
            numel(mitm_hand.FrameSubFileIndices{mitm_hand.Frame})>1);
        for i = 1:length(algorithm_info)
            algorithm_can_use_phd = ~strcmpi(algorithm_info(i).DataType,'Complex');
            algorithm_can_use_complex = ~strcmpi(algorithm_info(i).DataType,'Phase History');
            algorithm_must_have_multipol = strcmpi(algorithm_info(i).Polarimetric,'Yes');
            if strcmpi(algorithm_info(i).Selected,'Yes') && ...
                    ((data_is_phd && algorithm_can_use_phd) || ...
                    (data_is_complex && algorithm_can_use_complex && ...
                ~(~data_is_multipol && algorithm_must_have_multipol)) || ...
                    (~data_is_phd && ~data_is_complex)) % No data loaded, show all algorithms
                uimenu(algorithms_uimenu, 'Label', algorithm_info(i).TextName, ...
                    'Callback', @(obj, event) run_tool(algorithm_info(i)));
            end
        end
        uimenu(algorithms_uimenu, 'Label', 'Select Algorithms...', ...
            'Callback', @select_algorithms, 'Separator','on');
    end

    function select_algorithms(obj, event)
        algorithm_info = AlgorithmSelection(algorithm_info);
        populate_algorithm_menu();
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////