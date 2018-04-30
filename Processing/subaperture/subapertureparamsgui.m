function subapertureparamsgui(varargin)
%SUBAPERTUREPARAMSGUI Allows one to set subaperture parameters through a GUI
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
    'Name','Subaperture','Units','Pixels','Position',[0 0 300 250]);
movegui(fig_hand,'center');

uicontrol('Style','Text','String','Number of Frames:',...
    'Units','normalized','Position',[.1 .7 .2 .2],'parent',fig_hand,...
    'BackgroundColor',get(fig_hand,'Color'));
edit_frames = uicontrol('Style','Edit','String','7',...
    'Units','normalized','Position',[.3 .8 .15 .1],'parent',fig_hand);
uicontrol('Style','Text','String','Aperture Fraction:',...
    'Units','normalized','Position',[.1 .5 .2 .2],'parent',fig_hand,...
    'BackgroundColor',get(fig_hand,'Color'));
edit_ap = uicontrol('Style','Edit','String','0.25',...
    'Units','normalized','Position',[.3 .6 .15 .1],'parent',fig_hand);

bg_dir = uibuttongroup('Units','normalized','Position',[.1 .1 .35 .4],...
    'Title','Direction','BackgroundColor',get(fig_hand,'Color'));
radio_st = uicontrol('Style','Radio','String','Slow-Time',...
    'Units','normalized','Position',[.1 .65 .9 .2],'parent',bg_dir,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_ft = uicontrol('Style','Radio','String','Fast-Time',...
    'Units','normalized','Position',[.1 .25 .9 .2],'parent',bg_dir,...
    'BackgroundColor',get(fig_hand,'Color'));
bg_method = uibuttongroup('Units','normalized','Position',[.55 .35 .35 .55],...
    'Title','Method','BackgroundColor',get(fig_hand,'Color'));
radio_fp = uicontrol('Style','Radio','String','Full Pixel',...
    'Units','normalized','Position',[.2 .75 .8 .2],'parent',bg_method,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_norm = uicontrol('Style','Radio','String','Normal',...
    'Units','normalized','Position',[.2 .45 .8 .2],'parent',bg_method,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_min = uicontrol('Style','Radio','String','Minimal',...
    'Units','normalized','Position',[.2 .15 .8 .2],'parent',bg_method,...
    'BackgroundColor',get(fig_hand,'Color'));
uicontrol('Style','pushbutton','String','Run','FontWeight','bold',...
    'Units','normalized','Position',[.55 .1 .35 .15],'parent',fig_hand,...
    'Callback',@Run_Callback);

    %% Callback for button press
    function Run_Callback(hObject, eventdata)
        if get(radio_st,'Value')
            Dim = 1;
        else
            Dim = 2;
        end
        if get(radio_norm,'Value')
            method = 'normal';
        elseif get(radio_fp,'Value')
            method = 'fullpixel';
        else
            method = 'minimal';
        end
        
        if ~isempty(p.Results.filename)
            %now generate subapertures for specified file and AOI
            if (p.Results.aoi(1) == 0)
                AzLimits = 'full';
                RgLimits = 'full';
            else
                AzLimits = round([p.Results.aoi(1)  p.Results.aoi(1)+p.Results.aoi(3)-1]);
                RgLimits = round([p.Results.aoi(2)  p.Results.aoi(2)+p.Results.aoi(4)-1]);
            end
            tmp_name=tempname;
            subaperturefile(p.Results.filename,tmp_name,...
                'azlimits',AzLimits,'rnglimits',RgLimits,...
                'framenumber',p.Results.segment,...
                'apfraction',str2double(get(edit_ap,'String')),...
                'dim',Dim,...
                'frames',str2double(get(edit_frames,'String')),...
                'method',method);
            
            %launch MITM viewer to view
            outputfilestructs=dir([tmp_name '*']);
            for i=1:length(outputfilestructs)
                tmp_dir=fileparts(tmp_name);
                fulloutputfilenames{i}=[tmp_dir filesep outputfilestructs(i).name];
            end
            
            mitm_viewer(fulloutputfilenames);
            set(gcf,'DeleteFcn',@(obj,eventdata) file_cleanup(obj,eventdata,fulloutputfilenames)); % Delete temp files on closing
        else
            subaperturedemo('apfraction',str2double(get(edit_ap,'String')),...
                'dim',Dim,...
                'frames',str2double(get(edit_frames,'String')),...
                'method',method);
        end
        close(fig_hand);
    end
end

function file_cleanup(obj, eventdata, fulloutputfilenames)
    delete(get(obj,'Children')); % Calls destructor for hg_mitm_viewer, which releases files
    % Clean up intermediate files
    for j=1:length(fulloutputfilenames)
        delete(fulloutputfilenames{j});
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////