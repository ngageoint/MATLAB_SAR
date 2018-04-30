function imcontrast2(varargin)
% IMCONTRAST2 Simple replacement code for the IMCONTRAST tool in the image
% processing toolbox
%
% Author: Tim Cox, NRL
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Draw figure components
fig_hand = figure('Name','imcontrast2','NumberTitle','off',...
    'MenuBar','none','Units','characters','Position',[0 0 136 28]);
movegui(fig_hand,'center');

temp1 = uicontrol('Style','text','String','Image Cutoff Specification ','Parent',fig_hand,...
    'Units','normalized','Position',[0.3 0.9 0.4 0.06],...
    'FontSize',12,'FontWeight','bold','BackgroundColor',get(fig_hand,'Color'));
ax_hist = axes('Parent',fig_hand,...
    'Units','normalized', 'Position',[0.05 0.2 0.9 0.65]);
temp2 = uicontrol('Style','text','String','Min Value:','Parent',fig_hand,...
    'Units','characters','Position',[6 1.7 14.2 1.4],...
    'FontSize',10,'FontWeight','bold','BackgroundColor',get(fig_hand,'Color'))
edit_min = uicontrol('Style','edit','Parent',fig_hand,...
    'Units','characters','Position',[22.4 1.7 16.6 1.5]);
temp3 = uicontrol('Style','text','String','Max Value:','Parent',fig_hand,...
    'Units','characters','Position',[44.2 1.7 15.4 1.4],...
    'FontSize',10,'FontWeight','bold','BackgroundColor',get(fig_hand,'Color'));
edit_max = uicontrol('Style','edit','Parent',fig_hand,...
    'Units','characters','Position',[61 1.7 16.6 1.5]);
temp4 = uicontrol('Style','pushbutton','String','Update','Parent',fig_hand,...
    'Units','characters','Position',[92.6 1.5 17.4 1.7],...
    'Callback',@Update_Callback);

%% Setup histogram
axes_handle = varargin{1};
%set up min/max as current values
CLim = get(axes_handle,'CLim');
set(edit_min,'String',sprintf('%0.3g',CLim(1)));
set(edit_max,'String',sprintf('%0.3g',CLim(2)));

%plot histogram
cdata = get(findobj(axes_handle,'type','image'),'CData'); %get data
[n xout] = hist(double(cdata(:)),50);
bar(ax_hist,xout,n,1);
maxn = max(n);
set(ax_hist,'YLim',[0 maxn],'YTick',[]);

%draw min/max lines
hold(ax_hist,'on');
minh = plot(ax_hist,[CLim(1) CLim(1)],[0 maxn],'-r');
maxh = plot(ax_hist,[CLim(2) CLim(2)],[0 maxn],'-r');
hold(ax_hist,'off');


    %% Callback for update button
    function Update_Callback(hObject, eventdata)
        MinVal = str2double(get(edit_min,'String'));
        MaxVal = str2double(get(edit_max,'String'));
        
        %Make sure input is valid
        
        %draw min/max lines
        delete([minh maxh]);
        hold(ax_hist,'on');
        minh = plot(ax_hist,[MinVal MinVal],[0 maxn],'-r');
        maxh = plot(ax_hist,[MaxVal MaxVal],[0 maxn],'-r');
        hold(ax_hist,'off');
        
        %update CLim on original axis with new values
        set(axes_handle,'CLim',[MinVal MaxVal]);
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////