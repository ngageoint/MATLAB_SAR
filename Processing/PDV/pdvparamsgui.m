function pdvparamsgui(varargin)
%PDVPARAMSGUI Allows one to set PDV parameters through a GUI
%
% Written by: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input arguments
p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[0 0 0 0]);
p.addParamValue('segment',1);
p.parse(varargin{:});

%% Draw figure
fig_hand = figure('MenuBar','none','Toolbar','none','NumberTitle','off',...
    'Name','PDV','Units','Pixels','Position',[0 0 300 250]);
movegui(fig_hand,'center');

% Pixel shift
uicontrol('Style','Text','String','Pixel Shift:',...
    'Units','normalized','Position',[.05 .7 .2 .1],'parent',fig_hand,...
    'BackgroundColor',get(fig_hand,'Color'));
edit_ps = uicontrol('Style','Edit','String','0.25',...
    'Units','normalized','Position',[.25 .7 .15 .11],'parent',fig_hand);

% Direction
bg_dir = uibuttongroup('Units','normalized','Position',[.05 .25 .35 .35],...
    'Title','Direction','BackgroundColor',get(fig_hand,'Color'));
radio_st = uicontrol('Style','Radio','String','Slow-Time',...
    'Units','normalized','Position',[.1 .65 .9 .2],'parent',bg_dir,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_ft = uicontrol('Style','Radio','String','Fast-Time',...
    'Units','normalized','Position',[.1 .25 .9 .2],'parent',bg_dir,...
    'BackgroundColor',get(fig_hand,'Color'));

% Window size
uip_smooth = uipanel('Units','normalized','Position',[.45 .6 .5 .35],'parent',fig_hand,...
    'Title','Smoothing Window','BackgroundColor',get(fig_hand,'Color'));
uicontrol('Style','Text','String','Azimuth:',...
    'Units','normalized','Position',[.05 .6 .4 .3],'parent',uip_smooth,...
    'BackgroundColor',get(fig_hand,'Color'));
edit_az = uicontrol('Style','Edit','String','5',...
    'Units','normalized','Position',[.45 .6 .2 .35],'parent',uip_smooth);
uicontrol('Style','Text','String','Pixels',...
    'Units','normalized','Position',[.7 .6 .25 .3],'parent',uip_smooth,...
    'BackgroundColor',get(fig_hand,'Color'));
uicontrol('Style','Text','String','Range:',...
    'Units','normalized','Position',[.05 .1 .4 .3],'parent',uip_smooth,...
    'BackgroundColor',get(fig_hand,'Color'));
edit_rng = uicontrol('Style','Edit','String','5',...
    'Units','normalized','Position',[.45 .1 .2 .35],'parent',uip_smooth);
uicontrol('Style','Text','String','Pixels',...
    'Units','normalized','Position',[.7 .1 .25 .3],'parent',uip_smooth,...
    'BackgroundColor',get(fig_hand,'Color'));

% Fitler type
bg_ft = uibuttongroup('Units','normalized','Position',[.45 .25 .5 .35],...
    'Title','Filter Type','BackgroundColor',get(fig_hand,'Color'));
radio_mean = uicontrol('Style','Radio','String','Mean',...
    'Units','normalized','Position',[.1 .65 .9 .2],'parent',bg_ft,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_med = uicontrol('Style','Radio','String','Median',...
    'Units','normalized','Position',[.1 .25 .9 .2],'parent',bg_ft,...
    'BackgroundColor',get(fig_hand,'Color'));
uicontrol('Style','pushbutton','String','Run','FontWeight','bold',...
    'Units','normalized','Position',[.35 .05 .3 .15],'parent',fig_hand,...
    'Callback',@Run_Callback);

    %% Callback for button press
    function Run_Callback(hObject, eventdata)
        if get(radio_st,'Val')
            Dim = 1;
        else
            Dim = 2;
        end
        
        if get(radio_mean,'Value')
            filtertype = 'mean';
        else
            filtertype = 'median';
        end

        
        if ~isempty(p.Results.filename)
            %now generate PDV for specified file and AOI
            if (p.Results.aoi(1) == 0)
                AzLimits = 'full';
                RgLimits = 'full';
            else
                AzLimits = round([p.Results.aoi(1)  p.Results.aoi(1)+p.Results.aoi(3)-1]);
                RgLimits = round([p.Results.aoi(2)  p.Results.aoi(2)+p.Results.aoi(4)-1]);
            end
            tmp_name=tempname;
            tmp_name2=tempname;
            
            pdvfile(p.Results.filename,tmp_name,...
                'azlimits',AzLimits,'rnglimits',RgLimits,...
                'framenumber',p.Results.segment,...
                'deltax',str2double(get(edit_ps,'String')),...
                'dim',Dim,...
                'filtersize',[str2double(get(edit_az,'String')) str2double(get(edit_rng,'String'))],...
                'filtertype',filtertype);
            chipfile(p.Results.filename,tmp_name2,...
                'azlimits',AzLimits,'rnglimits',RgLimits,...
                'framenumber',p.Results.segment);
            
            %launch MITM viewer to view PDV
            mitm_viewer({tmp_name tmp_name2},'initialRemap','linearremap');
            set(gcf,'DeleteFcn',@file_cleanup); % Delete temp files on closing
        else
            pdvdemo('deltax',str2double(get(edit_ps,'String')),...
                'dim',Dim,...
                'filtersize',[str2double(get(edit_az,'String')) str2double(get(edit_rng,'String'))],...
                'filtertype',filtertype);
        end
        close(fig_hand);
        
        function file_cleanup(obj, eventdata)
            delete(get(obj,'Children')); % Calls destructor for hg_mitm_viewer, which releases files
            % Clean up intermediate files
            delete(tmp_name);
            delete(tmp_name2);
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////