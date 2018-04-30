function [ pixels_meta_agree ] = fs_vis_test( filename )
%FS_VIS_TEST Spatial frequency support visualization tool
%
% Allows a user to visualize potentially spatial variant spatial frequency
% support in SAR complex data and compare it against what is described in
% that file's metadata.  User can report back whether or not their manual
% review shows that the pixel data and metadata agree or not.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Variable setup
% Give output a default if window is closed with no buttons pressed
pixels_meta_agree = true;
% Constants for GUI layout
BORDER_SIZE = 20;
BUTTON_SIZE = [50 20];
BUFFER_SIZE = 50;
button_area_height = BUTTON_SIZE(2) + 3*BORDER_SIZE; % Area for instruction and buttons on bottom
AOI_SIZE = 512; % Generally a good idea to pick a power of 2, for faster FFTs
LINE_WIDTH = 2;  % We sometimes move this up for better PPT graphics
% Create "global" variables for overlays that subfunctions can see
row_DeltaKCOAPoly = 0; % Default until computed otherwise
col_DeltaKCOAPoly = 0;
row_bw_lines = [];
col_bw_lines = [];

%% Load and display an overview of the data
ro = open_reader(filename);
if iscell(ro), ro = ro{1}; end; % Don't account for multi-frame files
meta = ro.get_meta();
% Set defaults if not explicitly populated
if ~isfield(meta.Grid.Row,'DeltaKCOAPoly')
    meta.Grid.Row.DeltaKCOAPoly = 0;
end
if ~isfield(meta.Grid.Col,'DeltaKCOAPoly')
    meta.Grid.Col.DeltaKCOAPoly = 0;
end
% Two ways to pick bounds for our image axes sizes.  We assume two
% side-by-side images of equal size.
% 1) As large as will fit on screen:
screensize=get(0,'ScreenSize');
non_image_size = [BORDER_SIZE*5 BORDER_SIZE*2+button_area_height]; % Space we will need in figure beyond the two image axes
max_img_size = screensize(3:4)-screensize(1:2)+1;
max_img_size = max_img_size - non_image_size - BUFFER_SIZE; % Largest size that will fit on our screen
max_img_size(1) = floor(max_img_size(1)/2); % We need room for two side-by-side image panels
% 2) Arbitrary limit to keep small size for faster performance
% max_img_size = [512 512];
DecVal = max([1 ceil(double([meta.ImageData.NumCols meta.ImageData.NumRows])./max_img_size)]);
overview = densityremap(double(ro.read_chip([],[],[DecVal DecVal]))).';

%% Setup layout of figure and image axes
image_panel_size = min(max_img_size); % Make it square-- just because
% Figure
fig_size = [((2*image_panel_size)+(BORDER_SIZE*5)) ...
    (image_panel_size+(2*BORDER_SIZE)+button_area_height)];
fig_hand = figure('Name','Spatially variant frequency support', ...
    'NumberTitle','off','MenuBar','none','WindowStyle','modal', ...
    'Units','pixels','Position',[0 0 fig_size]);
movegui(fig_hand,'center');
% Axes
overview_axes_hand = axes('Parent',fig_hand,'Units','pixels', ...
    'Position',[BORDER_SIZE (BORDER_SIZE+button_area_height) ...
    image_panel_size image_panel_size]);
colormap(overview_axes_hand,gray(256));
transform_axes_hand = axes('Parent',fig_hand,'Units','pixels','NextPlot','add', ...
    'Position',[(image_panel_size+4*BORDER_SIZE) (BORDER_SIZE+button_area_height) ...
    image_panel_size image_panel_size]);
colormap(transform_axes_hand,gray(256));
image(overview,'Parent',overview_axes_hand,'CDataMapping','scaled');
set(overview_axes_hand, 'DataAspectRatio', [1 1 1], 'Visible', 'off');
transform_image_hand = image(0,... % Create empty MATLAB image object for later use
    'Parent',transform_axes_hand,'CDataMapping','scaled');
set(transform_axes_hand, 'Visible', 'on', ...
    'XLim',[0.5 AOI_SIZE+0.5], ...
    'YLim',[0.5 AOI_SIZE+0.5]);
% Axes titles
uicontrol('Style','text','String','Image Overview','Parent',fig_hand,...
    'Units','pixels','Position',[BORDER_SIZE ...
    (image_panel_size + button_area_height + BORDER_SIZE) image_panel_size BORDER_SIZE],...
    'FontSize',12,'FontWeight','bold','BackgroundColor',get(fig_hand,'Color'));
uicontrol('Style','text','String','Raw FFT of AOI','Parent',fig_hand,...
    'Units','pixels','Position',[(image_panel_size + (3*BORDER_SIZE)) ...
    (image_panel_size + button_area_height + BORDER_SIZE) image_panel_size BORDER_SIZE],...
    'FontSize',12,'FontWeight','bold','BackgroundColor',get(fig_hand,'Color'));

%% Setup AOI selector imrect
% Create imrect
init_start = floor(size(overview)/2 - AOI_SIZE/(2*DecVal)) + 1;
init_stop = floor(AOI_SIZE/DecVal);
init_start = max(init_start,1);
init_stop = min(init_stop, size(overview));
init_AOI = [init_start([2 1]) init_stop];
aoi_hand = imrect(overview_axes_hand,init_AOI);
api_handle = iptgetapi(aoi_hand);
fcn = makeConstrainToRectFcn('imrect',...
    get(overview_axes_hand,'XLim'),get(overview_axes_hand,'YLim'));
setPositionConstraintFcn(aoi_hand,fcn);
% Define callback function that will transform data into the spatial
% frequency domain and display.
if meta.Grid.Row.Sgn < 0 % Assume that Row.Sgn and Col.Sgn must be equal
    sf_fun = @(x) fftshift(fft2(x));
else
    sf_fun = @(x) ifftshift(ifft2(x));
end
api_handle.addNewPositionCallback(@(pos) newpos(sf_fun, pos));

%% Setup overlays of DeltaK1/DeltaK2
% We only need to do this once, not at each AOI update, since these are
% constant for the entire image.
krow_bounds = [round((AOI_SIZE-1) * ...
    (0.5 + meta.Grid.Row.SS*meta.Grid.Row.DeltaK1)) + 1, ...
    round((AOI_SIZE-1) * ...
    (0.5 + meta.Grid.Row.SS*meta.Grid.Row.DeltaK2)) + 1];
kcol_bounds = [round((AOI_SIZE-1) * ...
    (0.5 + meta.Grid.Col.SS*meta.Grid.Col.DeltaK1)) + 1, ...
    round((AOI_SIZE-1) * ...
    (0.5 + meta.Grid.Col.SS*meta.Grid.Col.DeltaK2)) + 1];
% For cases where DeltaKCOAPoly=0, these lines should be covered later by
% the DeltaKCOAPoly/ImpRespBW overlays.  However, for spatially variant
% frequency support, it may be important to see both bounds independently.
deltak_lines = plot(transform_axes_hand, ...
    xlim(transform_axes_hand), krow_bounds(1)*[1 1], '--r', ...
    xlim(transform_axes_hand), krow_bounds(2)*[1 1], '--r', ...
    kcol_bounds(1)*[1 1], ylim(transform_axes_hand), '--r', ...
    kcol_bounds(2)*[1 1], ylim(transform_axes_hand), '--r');
set(deltak_lines,'LineWidth',LINE_WIDTH);
newpos(sf_fun, init_AOI); % Initialize display by running once

%% Show instructions and buttons
% Allow user to confirm or deny agreement between pixel and XML metadata.
% Return confirmation as output argument
instructions = ['As AOI is moved across image, does the spatial frequency ' ...
    'support seen in the data agree with the bounds predicted from the ' ...
    'metadata (as shown by the graphics overlaid on the FFT image)?'];
uicontrol('Style','text','String',instructions,'Parent',fig_hand,...
    'Units','pixels','Position',[BORDER_SIZE BORDER_SIZE ...
    (2*image_panel_size-BUTTON_SIZE(1)-BORDER_SIZE) 2*BORDER_SIZE],...
    'FontSize',12,'HorizontalAlignment','left','BackgroundColor',get(fig_hand,'Color'));
uicontrol('Style','pushbutton','String','Yes','Parent',fig_hand,...
    'Units','pixels','Position',[(fig_size(1)-2*BORDER_SIZE-2*BUTTON_SIZE(1)) BORDER_SIZE BUTTON_SIZE],...
    'Callback',@yes_callback);
uicontrol('Style','pushbutton','String','No','Parent',fig_hand,...
    'Units','pixels','Position',[(fig_size(1)-BORDER_SIZE-BUTTON_SIZE(1)) BORDER_SIZE BUTTON_SIZE],...
    'Callback',@no_callback);
uiwait(fig_hand);
ro.close();

    function newpos(transform_fun, pos)
        % Display spatial frequency for new AOI
        % The size of data we read in is only dependent on the size of the
        % AOI_SIZE, and not on the size of the imrect in the
        % overview image.  We just use the center of the imrect, regardless
        % of its actual size.
        center_pos = (pos(1:2)+(pos(3:4)/2))*DecVal;
        start = max(floor(center_pos - (AOI_SIZE/2)),1);
        stop = min(start + AOI_SIZE, ...
            double([meta.ImageData.NumCols meta.ImageData.NumRows]));
        subimage = double(ro.read_chip([start(1) stop(1)],[start(2) stop(2)])).';
        set(transform_image_hand,'CData',logremap(transform_fun(subimage)));
        
        % Update overlay graphics for local DeltaKCOAPoly/ImpRespBW bounds
        if isfield(meta.Grid.Row, 'DeltaKCOAPoly') && ...
                any(meta.Grid.Row.DeltaKCOAPoly(:)) || isempty(row_bw_lines)
            row_DeltaKCOAPoly = sicd_polyval2d(meta.Grid.Row.DeltaKCOAPoly,...
                center_pos(1),center_pos(2),meta);
            row_bw = mod(round((AOI_SIZE-1) * ...
                (0.5 + meta.Grid.Row.SS * (row_DeltaKCOAPoly + ...
                ([-1 +1]*meta.Grid.Row.ImpRespBW/2)))) + 1,AOI_SIZE);
            delete(row_bw_lines); % Remove any previous markers
            row_bw_lines = plot(transform_axes_hand, ...
                xlim(transform_axes_hand), row_bw(1)*[1 1], '--b', ...
                xlim(transform_axes_hand), row_bw(2)*[1 1], '--b');
            set(row_bw_lines,'LineWidth',LINE_WIDTH);
            row_labels_all = {'-Kr_{BW}/2','Kr_{BW}/2','\DeltaKrow1','\DeltaKrow2','0'};
            row_vals_all = [row_bw krow_bounds mean(xlim(transform_axes_hand))];
            [row_vals, idx] = unique(row_vals_all);
            row_labels = row_labels_all(idx);
            set(transform_axes_hand, 'YTick', row_vals);
            set(transform_axes_hand, 'YTickLabel', row_labels);
        end
        if isfield(meta.Grid.Col, 'DeltaKCOAPoly') && ...
                any(meta.Grid.Col.DeltaKCOAPoly(:)) || isempty(col_bw_lines)
            col_DeltaKCOAPoly = sicd_polyval2d(meta.Grid.Col.DeltaKCOAPoly,...
                center_pos(1),center_pos(2),meta);
            col_bw = mod(round((AOI_SIZE-1) * ...
                (0.5 + meta.Grid.Col.SS * (col_DeltaKCOAPoly + ...
                ([-1 +1]*meta.Grid.Col.ImpRespBW/2)))) + 1,AOI_SIZE);
            delete(col_bw_lines); % Remove any previous markers
            col_bw_lines = plot(transform_axes_hand, ...
                col_bw(1)*[1 1], ylim(transform_axes_hand), '--b', ...
                col_bw(2)*[1 1], ylim(transform_axes_hand), '--b');
            set(col_bw_lines,'LineWidth',LINE_WIDTH);
            % Update axes ticks
            col_labels_all = {'-Kc_{BW}/2','Kc_{BW}/2','\DeltaKcol1','\DeltaKcol2','0'};
            col_vals_all = [col_bw kcol_bounds mean(xlim(transform_axes_hand))];
            [col_vals, idx] = unique(col_vals_all);
            col_labels = col_labels_all(idx);
            set(transform_axes_hand, 'XTick', col_vals);
            set(transform_axes_hand, 'XTickLabel', col_labels);
        end
    end

    function yes_callback(varargin)
        pixels_meta_agree = true;
        close(fig_hand);
    end

    function no_callback(varargin)
        pixels_meta_agree = false;
        close(fig_hand);
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////