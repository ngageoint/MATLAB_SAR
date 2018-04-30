function MIM( varargin )
%MIM MATLAB Icon Maker
%   GUI for creating SAR metadata icons
%
% Author: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

BACKGROUND_COLOR = [51 102 153]/255; % Default matches JIM

%% Setup GUI layout
h_figure = figure('Name','Matlab Icon Maker',... 
    'Position',[100 100 400 500],...
    'MenuBar','none',...
    'Units','normalized');
h_metaicon = axes('Parent',h_figure,...
    'Position',[.05 .25 .9 .7],...
    'Color',BACKGROUND_COLOR,...
    'XTick',[],...
    'YTick',[],...
    'Box','on');
BUTTON_HEIGHT = .08;
h_load_button = uicontrol('Parent',h_figure,...
    'Style','pushbutton',...
    'String','Load File',...
    'Units','normalized',...
    'Position',[.1 .14 .35 BUTTON_HEIGHT],...
    'Callback',@LoadFile_Callback);
h_save_button = uicontrol('Parent',h_figure,...
    'Style','pushbutton',...
    'String','Save Icon',...
    'Units','normalized',...
    'Position',[.55 .14 .35 BUTTON_HEIGHT],...
    'Enable','off',...
    'Callback',@SaveIcon_Callback);
h_batch_button = uicontrol('Parent',h_figure,...
    'Style','pushbutton',...
    'String','Batch Process',...
    'Units','normalized',...
    'Position',[.1 .03 .5 BUTTON_HEIGHT],...
    'Callback',@BatchProcess_Callback);
h_recursive = uicontrol('Parent',h_figure,...
    'Style','checkbox',...
    'String','Recursive',...
    'Units','normalized',...
    'Position',[.7 .03 .2 BUTTON_HEIGHT],...
    'BackgroundColor',get(h_figure,'Color'));

% Determine if data filename was passed as argument
p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.parse(varargin{:});

if length(p.Results.filename) > 0
    LoadImage(p.Results.filename);
end
   
%% Callbacks
    % --- Executes on button press in LoadFile.
    function LoadFile_Callback(hObject, eventdata)
        % hObject    handle to LoadFile (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        
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
        [fname, path] = uigetfile( sar_file_extensions({'complex','phd'}), 'Open SAR data file',pathstr,'MultiSelect', 'on' );
        if isnumeric(fname)
            return;
        end
        
        filename = strcat(path,fname);
        
        setpref('matlab_sar_toolbox','last_used_directory',path); %store path
        
        LoadImage(filename);
    end

    % Load image helper function
    function LoadImage(filename)
        % determine if input file is complex image or phase history data
        if ~isempty(guess_ph_format(filename))
            MetaIcon_PHD(filename,'handle',h_metaicon,'allow_editing',true,...
                'background_color', BACKGROUND_COLOR);
        else
            MetaIcon_Complex(filename,'handle',h_metaicon,'allow_editing',true,...
                'background_color', BACKGROUND_COLOR);
        end
        
        set(h_save_button,'Enable','on');
    end

    % --- Executes on button press in SaveIcon.
    function SaveIcon_Callback(hObject, eventdata)
        % hObject    handle to SaveIcon (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        
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

        formats = imformats; % All image formats handled by MATLAB
        formats = formats(~cellfun(@isempty,{formats.write})); % All writable formats
        extensions =  cellfun(@(x) x{1},{formats.ext},'UniformOutput',false)';
        filters = cellfun(@(x) ['*.' x],extensions,'UniformOutput',false);
        [fname, path, FilterIndex] = uiputfile( [filters, {formats.description}'],'Save File',pathstr);
        
        if fname
            filename = strcat(path,fname);
            
            setpref('matlab_sar_toolbox','last_used_directory',path); %store path
            
            frame = getframe(h_metaicon);
            
            imwrite(frame.cdata,filename,extensions{FilterIndex});
        end
    end

    % --- Executes on button press in BatchProcess.
    function BatchProcess_Callback(hObject, eventdata)
        % hObject    handle to BatchProcess (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        
        %allow user to select directory
        %load last path
        if ispref('matlab_sar_toolbox','last_used_directory')
            pathstr = getpref('matlab_sar_toolbox','last_used_directory');
            if ~ischar(pathstr)||~exist(pathstr,'dir')
                pathstr = pwd;
            end
        else
            pathstr = pwd;
        end

        infolder = uigetdir(pathstr,'Select Input Directory');
        outfolder = uigetdir(infolder,'Select Output Directory');
        
        setpref('matlab_sar_toolbox','last_used_directory',infolder); %store path
        
        %search for all image & phd files (we may want to do support files at a later
        %time)
        cell_for_uigetfile = sar_file_extensions('complex');
        ext_filters = textscan(cell_for_uigetfile{1,1}','%s',Inf,'delimiter',';');
        complex_extensions = cellfun(@(x) x(2:end), ext_filters{1},'UniformOutput',false); % Remove '*'
        cell_for_uigetfile = sar_file_extensions('phd');
        ext_filters = textscan(cell_for_uigetfile{1,1}','%s',Inf,'delimiter',';');
        phd_extensions = cellfun(@(x) x(2:end), ext_filters{1},'UniformOutput',false); % Remove '*'
        extensions = [complex_extensions; phd_extensions] ;
        
        if get(h_recursive,'Value')
            FileNames = rdir(infolder,extensions);
            for i=1:length(FileNames)
                filename = FileNames{i};
                try
                    if any(strcmpi(phd_extensions, extensions{i}))
                        MetaIcon_PHD(filename,'handle',h_metaicon,...
                            'background_color', BACKGROUND_COLOR);
                    else
                        MetaIcon_Complex(filename,'handle',h_metaicon,...
                            'background_color', BACKGROUND_COLOR);
                    end
                    %if icon generation was sucessful, then save png
                    [path, fname] = fileparts(filename);
                    frame = getframe(h_metaicon);
                    imwrite(frame.cdata,fullfile(outfolder,[fname '.png']),'png');
                catch
                end
            end
        else
            for i=1:length(extensions)
                %get list of files for current extension
                dirsearch = sprintf('%s/*%s',infolder,extensions{i});
                listing = dir(dirsearch);
                for j=1:length(listing)
                    filename = sprintf('%s/%s',infolder,listing(j).name);
                    try
                        if any(strcmpi(phd_extensions, extensions{i}))
                            MetaIcon_PHD(filename,'handle',h_metaicon,...
                                'background_color', BACKGROUND_COLOR);
                        else
                            MetaIcon_Complex(filename,'handle',h_metaicon,...
                                'background_color', BACKGROUND_COLOR);
                        end
                        %if icon generation was sucessful, then save png
                        [path, fname] = fileparts(listing(j).name);
                        frame = getframe(h_metaicon);
                        imwrite(frame.cdata,fullfile(outfolder,[fname '.png']),'png');
                    catch
                    end
                end
            end
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////