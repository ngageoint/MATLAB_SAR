function csiparamsgui(varargin)
%CSIPARAMSGUI Allows one to set CSI parameters through a GUI
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
    'Name','CSI','Units','Pixels','Position',[0 0 300 200]);
movegui(fig_hand,'center');

bg_dir = uibuttongroup('Units','normalized','Position',[.05 .5 .35 .4],...
    'Title','Direction','BackgroundColor',get(fig_hand,'Color'));
radio_st = uicontrol('Style','Radio','String','Slow-Time',...
    'Units','normalized','Position',[.1 .65 .9 .2],'parent',bg_dir,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_ft = uicontrol('Style','Radio','String','Fast-Time',...
    'Units','normalized','Position',[.1 .25 .9 .2],'parent',bg_dir,...
    'BackgroundColor',get(fig_hand,'Color'));
bg_legend = uibuttongroup('Units','normalized','Position',[.45 .35 .5 .55],...
    'Title','CSI Legend','BackgroundColor',get(fig_hand,'Color'));
check_leg = uicontrol('Style','checkbox','String','Create CSI Legend',...
    'Units','normalized','Position',[.1 .7 .9 .2],'parent',bg_legend,...
    'BackgroundColor',get(fig_hand,'Color'),'Value',false);
radio_nu = uicontrol('Style','Radio','String','North Up',...
    'Units','normalized','Position',[.2 .45 .8 .2],'parent',bg_legend,...
    'BackgroundColor',get(fig_hand,'Color'));
radio_ed = uicontrol('Style','Radio','String','Energy Down',...
    'Units','normalized','Position',[.2 .15 .8 .2],'parent',bg_legend,...
    'BackgroundColor',get(fig_hand,'Color'));
check_full_res = uicontrol('Style','checkbox','String','Full Resolution',...
    'Units','normalized','Position',[.1 .3 .35 .2],'parent',fig_hand,...
    'BackgroundColor',get(fig_hand,'Color'),'Value',true);
uicontrol('Style','pushbutton','String','Run','FontWeight','bold',...
    'Units','normalized','Position',[.35 .1 .3 .15],'parent',fig_hand,...
    'Callback',@Run_Callback);

    %% Callback for button press
    function Run_Callback(hObject, eventdata)
        if get(radio_st,'Val')
            Dim = 1;
        else
            Dim = 2;
        end
        
        if ~isempty(p.Results.filename)
            %now generate CSI for specified file and AOI
            if (p.Results.aoi(1) == 0)
                AzLimits = 'full';
                RgLimits = 'full';
            else
                AzLimits = round([p.Results.aoi(1)  p.Results.aoi(1)+p.Results.aoi(3)-1]);
                RgLimits = round([p.Results.aoi(2)  p.Results.aoi(2)+p.Results.aoi(4)-1]);
            end
            tmp_name=tempname;
            csifile(p.Results.filename,tmp_name,'azlimits',AzLimits,'rnglimits',RgLimits,...
                'Dim',Dim,'fullres',get(check_full_res,'Value'),'framenumber',p.Results.segment);
            
            %CSILegend
            if get(check_leg,'Value')
                if Dim ==1
                    if (p.Results.aoi(1) == 0)
                        args = {}
                    else
                        args = {[p.Results.aoi(1)+round(p.Results.aoi(3)/2),...
                            p.Results.aoi(2)+round(p.Results.aoi(4)/2)]};
                    end
                    CSILegend_Slow(p.Results.filename, get(radio_ed,'Value'), args{:});
                else
                    CSILegend_Fast(p.Results.filename);
                end
            end
            
            %launch MITM viewer to view CSI
            outputfilestructs=dir([tmp_name '*']);
            for i=1:length(outputfilestructs)
                tmp_dir=fileparts(tmp_name);
                fulloutputfilenames{i}=[tmp_dir filesep outputfilestructs(i).name];
            end
            
            mitm_viewer([tmp_name '.mbw']);
            set(gcf,'DeleteFcn',@(obj,eventdata) file_cleanup(obj,eventdata,fulloutputfilenames)); % Delete temp files on closing
        else
            csidemo('Dim',Dim,'fullres',get(check_full_res,'Value'));
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