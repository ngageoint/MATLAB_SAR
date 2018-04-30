function DrawMetaIcon(meta, varargin)
%DRAWMETAICON Draw an icon containing information from a SAR collect 
%
% DRAWMETAICON (meta, 'PropertyName', PropertyValue, ...) generates a
% standard icon representing a wide variety of SAR collects.
%
% INPUTS:
%   meta      - required: metadata structure with the following fields:
%      IID: Unique identified for image
%      Pol: 2 character description of polarization
%      ZuluTime: Start time (zulu) of collect in MATLAB serial date number
%      CC: (optional) Country code
%      IPR: Impulse response width (meters)
%      Duration: length of the collect (seconds)
%      EffectiveDuration: (optional) length of collection for a given point
%         in the scene (seconds)
%      SCP: Scene center point (meters)
%      ARP: Aperture reference point
%      ARV: Aperture velocity
%   Property name       Description
%      allow_editing    Whether or not to allow manual editing of text
%                       fields.  (Default is false.)
%      background_color Color of the background (in RGB).
%      content_color    Color of arrows, text, and icon (in RGB).
%      handle           Handle to MATLAB axes where the metaicon will be drawn
%      linewidth        The thickness of the line used to draw the arrows
%
% OUTPUTS:
%   none - icon is drawn to specified axis 
%
% Authors: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('allow_editing', false);
p.addParamValue('background_color', [51 102 153]/255); % Default matches JIM
p.addParamValue('content_color', [1 1 1]); % Default is white.
p.addParamValue('handle', gca());
p.addParamValue('linewidth', 2);
p.FunctionName = mfilename;
p.parse(varargin{:});
% For simpler notation
handle = p.Results.handle;
linewidth = p.Results.linewidth;
content_color = p.Results.content_color;

%% Setup axes
cla(handle);
axis(handle, 'image');
axis(handle, 'ij');
axis(handle, [0, 1, 0, 1]);
set(handle, 'XTick', []);
set(handle, 'YTick', []);
set(handle, 'YDir', 'normal');
box(handle, 'on');
hold(handle, 'on');
set(handle, 'Color',p.Results.background_color);
% set(handle, 'Units', 'normalized');
% set(handle, 'Position', position);
% Plot to bounds so matlab doesn't get confused and try to re-scale (These
% lines won't be visible though).
plot(handle, [0 1 1 0], [0 0 1 1], 'k', 'LineWidth', 1);
    
%% Derive geometry info
g = vect2geom(meta.SCP,meta.ARP,meta.ARV);
LatLon = ecf_to_geodetic(meta.SCP);

%% Draw arrows
% Right/left arrow
plot(handle, [0.1 0.9], [0.05 0.05], 'Color', content_color, 'LineWidth', linewidth);
if g.right >= 0
    plot(handle, [0.9, 0.87], [.05, .07], 'Color', content_color, 'LineWidth', linewidth);
    plot(handle, [0.9, 0.87], [.05, .03], 'Color', content_color, 'LineWidth', linewidth);
else
    plot(handle, [0.1, 0.13], [.05, .07], 'Color', content_color, 'LineWidth', linewidth);
    plot(handle, [0.1, 0.13], [.05, .03], 'Color', content_color, 'LineWidth', linewidth);
end
 
% Ascend/descend arrow
sensor_altitude = norm(meta.ARP) - 6371000; % Distance from radius of earth
spaceborne = sensor_altitude > 80000; % 50 mile threshold for airborne vs. spaceborne
if spaceborne % Ascending/descending is only relevant to orbits
    plot(handle, [0.95, 0.95], [0.1, 0.9], 'Color', content_color, 'LineWidth', linewidth);
    if g.ascend >= 0
        plot(handle, [0.95, 0.93], [0.9, 0.87], 'Color', content_color, 'LineWidth', linewidth);
        plot(handle, [0.95, 0.97], [0.9, 0.87], 'Color', content_color, 'LineWidth', linewidth);
    else
        plot(handle, [0.95, 0.93], [0.1, 0.13], 'Color', content_color, 'LineWidth', linewidth);
        plot(handle, [0.95, 0.97], [0.1, 0.13], 'Color', content_color, 'LineWidth', linewidth);
    end
end

%% Show icon
xpos = [.6 .9]; ypos = [.5 .2]; % Where to place image
path = fileparts(mfilename('fullpath'));
if spaceborne
    sensor_mask = imread(fullfile(path, 'space.png'));
else
    sensor_mask = imread(fullfile(path, 'air.png'));
    xpos = xpos + 0.05; % Don't need to make room for ascend/descend arrow
end
if g.right < 0
    sensor_mask = fliplr(sensor_mask);
end
% Data read from file will be used as for transparency, not used as image
sensor_mask = double(sensor_mask);
transparency = repmat(sensor_mask/max(sensor_mask(:)),[1 1 3]);
content = repmat(permute(content_color(:),[3 2 1]),size(sensor_mask));
background = repmat(permute(p.Results.background_color(:),[3 2 1]),size(sensor_mask));
icon_with_transparency = content.*transparency + background.*(1-transparency);
h_image = imagesc(icon_with_transparency,'Parent',handle);
set(h_image, 'Xdata', xpos);
set(h_image, 'Ydata', ypos);
% The following line would have been the easier way to do it.
% Unfortunately, using transparency in an axes makes it unable to be
% extracted by getframe later, so we have to make or own transpency...
% set(h_image, 'AlphaData', sensor_mask/max(sensor_mask(:)));

%% Insert text
set(handle, 'DefaulttextColor', content_color);
LEFT_POS = 0.1;
% IID
h_iid = text(LEFT_POS, 0.85, sprintf('%s %s', meta.IID, meta.Pol), ...
    'Interpreter', 'none', 'Parent', handle);
% Time
try
    ZuluStr = [datestr(meta.ZuluTime,'HH:MM') ' Z']; % Zulu time in string format
    days_past_zulu=LatLon(2)/360;
    LocalStr = [datestr(meta.ZuluTime+days_past_zulu,'HH:MM') ' LS']; % Local solar time
catch
    ZuluStr = '';
    LocalStr = '';
end
h_time = text(LEFT_POS, 0.75, sprintf('Time: %s / %s', ZuluStr, LocalStr),...
    'Parent', handle);
% Geographic info
LatStr = latlonstr(LatLon(1),'lat','num_units',2);
LonStr = latlonstr(LatLon(2),'lon','num_units',2);
geo_str = sprintf('GEO: %s %s', LatStr, LonStr);
if isfield(meta,'CC')
    geo_str = sprintf('%s / %s', geo_str, meta.CC);
end
h_geo = text(LEFT_POS, 0.65, geo_str, 'Parent', handle);
% IPR
if isfield(meta,'units')&&strcmpi(meta.units,'metric')
    ipr_val = meta.IPR;
    ipr_units = 'm';
    if ipr_val<1
        ipr_val = ipr_val*100;
        ipr_units = 'cm';
    end
else
    ipr_val = meta.IPR/FEET_TO_METERS;
    ipr_units = 'ft';
    if ipr_val<1
        ipr_val = ipr_val*12;
        ipr_units = 'in';
    end
end
ipr_str = sprintf('IPR: %.1f %s / %.1f s', ipr_val, ipr_units, meta.Duration);
if isfield(meta,'EffectiveDuration')
    ipr_str = sprintf('%s / %.1f s', ipr_str, meta.EffectiveDuration);
end
h_ipr = text(LEFT_POS, 0.55, ipr_str, 'Parent', handle);
% Geometry info
degree_char = char(176);
h_graze = text(LEFT_POS, 0.45, sprintf('Graze: %+.2f%c', g.graze*180/pi, degree_char), ...
    'Parent',handle);
h_az = text(LEFT_POS, 0.35, sprintf('Azimuth: %+.2f%c', g.azimuth*180/pi, degree_char),...
    'Parent',handle);
h_layover = text(LEFT_POS, 0.25, sprintf('Layover: %+.2f%c', g.layover*180/pi, degree_char),...
    'Parent',handle);
h_multi = text(LEFT_POS, 0.15, sprintf('Multipath: %+.2f%c', g.multipath*180/pi, degree_char),...
    'Parent',handle);

%% Scale text to fit properly
all_text = [h_iid h_time h_geo h_ipr h_graze h_az h_layover h_multi];
set(all_text, 'VerticalAlignment', 'baseline');
set(all_text, 'FontUnits', 'normalized'); % Allow font size to change with axes resize
if p.Results.allow_editing
    set(all_text, 'ButtonDownFcn', @(src,evt) set(src,'Editing','on')); % Allow user to edit fields after drawn
    set(handle, 'ButtonDownFcn', @(src,evt) fit_text(all_text)); % Update text size after user edit
end
fit_text(all_text);

end

function fit_text(text_handles)
    DEFAULT_HEIGHT = 0.09; % Assures text fits in axes vertically
    MAX_WIDTH = 0.75;
    set(text_handles, 'FontSize', DEFAULT_HEIGHT);
    extent = cell2mat(get(text_handles, 'Extent')); % Extent of each text portion
    widest = max(extent(:,3));
    if widest>MAX_WIDTH % Longest line collides with ascend/descend arrow
        set(text_handles, 'FontSize', DEFAULT_HEIGHT * (MAX_WIDTH / widest)); % Assume all text the same size
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////