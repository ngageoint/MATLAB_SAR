function taser_about( )
%TASER_ABOUT Brings up the "About" window in the TASER menu
%
% Author: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fig_hand=figure('Name', 'Taser: About', 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'Toolbar', 'none');

% Draw title
text_ax = axes('Parent',fig_hand,'Position',[0.05 .7 0.5 .3]);
axis(text_ax, 'off');
title_text = text(.5, 0.6, 'TASER', 'FontSize', 20, 'FontWeight', 'bold', 'FontAngle', 'italic', 'HorizontalAlignment', 'center', 'Parent', text_ax);
if isdeployed()&&exist('taser_compile_date','file')
    set(title_text,'Position',[.5 .7]); % Make room for version text
    text(.5, 0.4, ['Version ' taser_compile_date], 'FontSize', 12, 'HorizontalAlignment', 'center', 'Parent', text_ax)
end
text(.5, 0.2, 'Naval Research Laboratory', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Parent', text_ax)
text(.5, 0.05, 'National Geospatial-Intelligence Agency', 'FontSize', 10, 'HorizontalAlignment', 'center', 'Parent', text_ax)

% Draw scroll bar with text info
license_ui = uicontrol(fig_hand, 'Style', 'edit',...
    'Units', 'normalized', ...
    'Position', [.05 .05 .55 .6], ...
    'min',0,'max',5,... % No effect, except max must be greater than min
    'HorizontalAlign','left',...
    'BackgroundColor',[1 1 1],...
    'enable','inactive');
scrollbar_string = '';
credits_path = fullfile(fileparts(mfilename('fullpath')),'taser_credits.txt');
if exist(credits_path,'file')
    fid = fopen(credits_path);
    scrollbar_string = fread(fid,Inf,'*char').';
    fclose(fid);
end
% Assumes this file is in Tools/MatlabWorkbench/about path of MATLAB SAR Toolbox
license_path = fullfile([fileparts(mfilename('fullpath')) filesep '..' filesep '..' filesep '..'],'license.txt');
if exist(license_path,'file')
    fid = fopen(license_path);
    fgetl(fid);
    license_strings = fread(fid,Inf,'*char').';
    fclose(fid);
    % Concatenate so that uicontrol will wrap properly.
    license_strings = regexprep(license_strings, '\r\n(?!\r\n)(?<!\r\n\r\n)', ''); % Single occurence only
    scrollbar_string = [scrollbar_string sprintf('\r\n\r\nLicense:\r\n\r\n') license_strings];
end
set(license_ui,'String',scrollbar_string);
             
% Draw seals
nga_ax = axes('Parent', fig_hand,'Position', [.6 .55 .4 .4]);
nrl_ax = axes('Parent', fig_hand,'Position', [.6 .05 .4 .4]);
image(imread(fullfile(fileparts(mfilename('fullpath')),'NGA_Seal.png'),...
    'BackgroundColor',get(fig_hand,'Color')),'Parent',nga_ax);
image(imread(fullfile(fileparts(mfilename('fullpath')),'NRL_Seal.png'),...
    'BackgroundColor',get(fig_hand,'Color')),'Parent',nrl_ax);
set([nga_ax nrl_ax], 'DataAspectRatio', [1 1 1]);
axis([nga_ax nrl_ax], 'off');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////