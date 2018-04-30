function f = MetaIcon(filename,varargin)
%METAICON - Function to draw a metaicon 
%
% Function to draw a meta-icon.  The metaicon is populated adaptively based 
% on the available data.  A filename or SICD meta structure can be passed
% in as the first argument.  This function draws vectors for layover,
% multipath, shadow and north directions.  The vectors are projected either
% in the slant plane (default) or optionally in the ground plane.
%
% INPUTS:
%   filename      - optional : filename or complex image meta structure (default=file browser) 
%   handle        - optional : handle to draw metaicon, (default=new figure)
%   GroundProject - optional : overide to plot vectors in ground plane (if
%                              image is in ground it will plot as such)
%   segment       - optional : image number of multi-image file (filename must be specified)
%
% OUTPUTS:
%   figure handle if new figure is generated, otherwise f is empty
%
% VERSION:
%   1.0
%     - Tim Cox 20170606
%     - initial version, based on previous meta icon functions written by
%       Tim Cox (NRL),Wade Schwartzkopf (NGA-R), Ralph Fiedler (NRL) and
%       Bob Jansen(NRL) 
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('handle',[]);
p.addParamValue('GroundProject',false);
p.addParamValue('segment',1);
p.parse(varargin{:});

%% Get input data
% If there is no filename/meta passed in then open a file browser and let
% the user select the filename
if ~exist('filename','var')
    % Load last path
    if ispref('matlab_sar_toolbox','last_used_directory')
        pathstr = getpref('matlab_sar_toolbox','last_used_directory');
        if ~ischar(pathstr)||~exist(pathstr,'dir')
            pathstr = pwd;
        end
    else
        pathstr = pwd;
    end
    
    [fname,pathstr]=uigetfile(sar_file_extensions({'complex','phd'}), ...
                        'Open Data File',pathstr,'MultiSelect','off');
                    
    filename = [pathstr fname];
    if ~fname, return; end  % Cancel
    setpref('matlab_sar_toolbox','last_used_directory',pathstr);
end

% Determine input data type
if isstruct(filename)
    meta = filename;
else
    if ~isempty(guess_ph_format(filename))
        reader_obj = open_ph_reader(filename);
        meta = reader_obj.get_meta();
        if ~isfield(meta,'SCPCOA')
            % Get info for a few specific pulses
            if isfield(reader_obj,'read_raw')
                read_fun = reader_obj.read_raw;
            elseif isfield(reader_obj,'read_cphd')
                read_fun = reader_obj.read_cphd;
            end
            [~, vbmeta] = read_fun([1 ...
                round(meta.Data.ArraySize(1).NumVectors/2) ...
                meta.Data.ArraySize(1).NumVectors],[]);
        end
        reader_obj.close();
    elseif ~isempty(guess_complex_format(filename))
        reader_obj = open_reader(filename);
        if iscell(reader_obj); reader_obj = reader_obj{p.Results.segment}; end
        meta = reader_obj.get_meta();
        reader_obj.close();
    end
end

%% Metadata Entries
% Collect timing
if isfield(meta,'Timeline')
    if isfield(meta.Timeline,'CollectStart')
        CollectStart = meta.Timeline.CollectStart;
    end
    if isfield(meta.Timeline,'CollectDuration')
        CollectDuration = meta.Timeline.CollectDuration;
    end
elseif isfield(meta,'Global')
    if isfield(meta.Global,'CollectStart')
        CollectStart = meta.Global.CollectStart;
    end
    if isfield(meta.Global,'CollectDuration')
        CollectDuration = meta.Global.CollectDuration;
    end
elseif exist('vbmeta','var')
    CollectDuration = vbmeta.TxTime(end);
end

% IID
if exist('CollectStart','var') && isfield(meta,'CollectionInfo') && ...
        isfield(meta.CollectionInfo,'CollectorName')
    %Use Collect Date and Collector Name if available
    IID_line = {[upper(datestr(CollectStart,'ddmmmyy')) ' ' ...
        meta.CollectionInfo.CollectorName(1:min(4,end))]};
elseif isfield(meta,'CollectionInfo') && isfield(meta.Global,'CoreName')
    IID_line = {meta.CollectionInfo.CoreName(1:min(16,end))};
else
    IID_line = {};
end

if isfield(meta,'GeoData') && isfield(meta.GeoData,'SCP') && isfield(meta.GeoData.SCP,'ECF')
    SCP = [meta.GeoData.SCP.ECF.X; meta.GeoData.SCP.ECF.Y; meta.GeoData.SCP.ECF.Z];
elseif isfield(meta,'SRP') 
    if strcmp(meta.SRP.SRPType,'FIXEDPT')
        SCP = meta.SRP.FIXEDPT.SRPPT;
        SCP = [SCP.X; SCP.Y; SCP.Z];
    elseif strcmp(meta.SRP.SRPType,'PVTPOLY')
        SCP = [polyval(meta.SRP.PVTPOLY.SRPPVTPoly.X(end:-1:1), meta.Global.CollectDuration);...
            polyval(meta.SRP.PVTPOLY.SRPPVTPoly.Y(end:-1:1), meta.Global.CollectDuration);...
            polyval(meta.SRP.PVTPOLY.SRPPVTPoly.Z(end:-1:1), meta.Global.CollectDuration)];
    end
elseif exist('vbmeta','var')
    SCP = vbmeta.SRPPos(2,:).';
end
if exist('SCP','var')
    lla = ecf_to_geodetic(SCP);
    Lat = lla(1);
    Lon = lla(2);
end

if exist('CollectStart','var')
    IID_line = [IID_line datestr(CollectStart,'HHMMZ')];
    if exist('SCP','var')
        try
            [LocalTime,Offset] = GetLocalTime(Lat,Lon,CollectStart);
            % LocalStr = [datestr(LocalTime,'HHMM'),'L']; % Local Time
            IID_line = [IID_line [num2str(Offset,'%+g') 'L']];
        end
    end
end 

try % CountryCode field or ConvertCountryCode function might not exist
    CC = meta.CollectionInfo.CountryCode(1:2);
    if iscell(CC), CC = CC{1}; end % Use the primary country if multiple
    CC = ConvertCountryCode(CC,'Trigraph');
catch %if no CC or lat/lon, try to get country from lat/lon
    try % GetCountryCodeName function might not exist of Lat/Lon might be empty
        CC = GetCountryCodeName([Lat Lon],'Trigraph');
    end
end

% Determine if Tgt info is in GeoData->GeoInfo
if isfield(meta,'GeoData') && isfield(meta.GeoData,'GeoInfo')
    if ~iscell(meta.GeoData.GeoInfo) && isscalar(meta.GeoData.GeoInfo)
        meta.GeoData.GeoInfo = {meta.GeoData.GeoInfo};
    end
    pri_index = cellfun(@(x) isfield(x,'name')&&strcmp(x.name,'PRIMARY_TARGET'),meta.GeoData.GeoInfo);
    if sum(pri_index)==1 && ...
            isfield(meta.GeoData.GeoInfo{pri_index}, 'Desc') && ...
            isfield(meta.GeoData.GeoInfo{pri_index}.Desc,'value')
        Tgt = meta.GeoData.GeoInfo{pri_index}.Desc.value(1:min(10,end));
    end
end

if isfield(meta,'ImageFormation') && isfield(meta.ImageFormation,'TxRcvPolarizationProc')
    Pol = meta.ImageFormation.TxRcvPolarizationProc; % Probably string, but could be cell array
elseif isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'RcvChannels') && ...
       isfield(meta.RadarCollection.RcvChannels,'ChanParameters') && ...
       isfield(meta.RadarCollection.RcvChannels.ChanParameters,'TxRcvPolarization')
    Pol = {meta.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization}; % RcvChannels is struct array
end
if exist('Pol','var')
    if iscell(Pol)
        Pol = cellfun(@(x) x(1:2:3), Pol, 'UniformOutput', false);
        if numel(Pol)==4
            Pol = 'QUAD';
        elseif numel(Pol)==2 && Pol{1}(1) == 'H' && Pol{2}(1) == 'H'
            Pol = 'HD';
        elseif numel(Pol)==2 && Pol{1}(1) == 'V' && Pol{2}(1) == 'V'
            Pol = 'VD';
        elseif numel(Pol)==1
            Pol = Pol{1};
        else % Don't know what to do.  Could be HH/VV or other unusual pol
            clear Pol;
        end
    else
        Pol = Pol(1:2:3);
    end
end

% See if we can find the RNIIRS
if isfield(meta,'CollectionInfo') && isfield(meta.CollectionInfo,'Parameter')
    if ~iscell(meta.CollectionInfo.Parameter) && isscalar(meta.CollectionInfo.Parameter)
        meta.CollectionInfo.Parameter = {meta.CollectionInfo.Parameter};
    end
    rniirs_index = cellfun(@(x) isfield(x,'name')&&strcmp(x.name,'PREDICTED_RNIIRS'), ...
        meta.CollectionInfo.Parameter);
    if sum(rniirs_index)==1 && ...
            isfield(meta.CollectionInfo.Parameter{rniirs_index}, 'value')
        RNIIRS = meta.CollectionInfo.Parameter{rniirs_index}.value;
    end
end

if isfield(meta,'SCPCOA')
    ARP = [meta.SCPCOA.ARPPos.X; meta.SCPCOA.ARPPos.Y; meta.SCPCOA.ARPPos.Z];
    ARV = [meta.SCPCOA.ARPVel.X; meta.SCPCOA.ARPVel.Y; meta.SCPCOA.ARPVel.Z];
elseif exist('vbmeta','var') && exist('SCP','var')
    % Sometimes no SCPCOA (non-standard field) for phase history
    ARP = (vbmeta.TxPos(2,:) + vbmeta.RcvPos(2,:))/2;
    ARV = ((((vbmeta.TxPos(end,:) + vbmeta.RcvPos(end,:))/2) - ...
        ((vbmeta.TxPos(1,:) + vbmeta.RcvPos(1,:))/2)) / ...
        (((vbmeta.TxTime(end,:) + vbmeta.RcvTime(end,:))/2) - ...
        ((vbmeta.TxTime(1,:) + vbmeta.RcvTime(1,:))/2))).';
    meta.GeoData.SCP.ECF = struct('X',SCP(1),'Y',SCP(2),'Z',SCP(3));
    meta.SCPCOA.ARPPos = struct('X',ARP(1),'Y',ARP(2),'Z',ARP(3));
    meta.SCPCOA.ARPVel = struct('X',ARV(1),'Y',ARV(2),'Z',ARV(3));
    meta = derived_sicd_fields(meta); % Compute SCPCOA angles.  Could also use vect2geom here, but this uses structure already in SICD
end

try % Lots of fields required for this that might not be there
    Lambda = 2/meta.Grid.Row.KCtr;
    Theta = meta.Grid.Col.ImpRespBW*Lambda/2;
    R = norm(abs(SCP-ARP));    
    V = norm(ARV);
    rov = (R/V);
    eff_ap_time = abs((Theta*rov)/sind(meta.SCPCOA.DopplerConeAng));
end

if isfield(meta,'Grid') %maintian resolution of plane that image is in...
    AzIPR = meta.Grid.Col.ImpRespWid/0.3048; % Meters to feet
    RgIPR = meta.Grid.Row.ImpRespWid/0.3048;     
    if p.Results.GroundProject && strcmp(meta.Grid.ImagePlane,'SLANT')
        AzIPR = AzIPR/cosd(meta.SCPCOA.TwistAng);
        RgIPR = RgIPR/cosd(meta.SCPCOA.GrazeAng);
    end
else
    if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
            isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
            isfield(meta.RadarCollection.Waveform.WFParameters,'TxRFBandwidth')
        BW = meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth/1e6; %MHz
    end
end

if isfield(meta,'SCPCOA')
    Azimuth = meta.SCPCOA.AzimAng;
    Layover = meta.SCPCOA.LayoverAng;
    MultipathGround = -atan(tand(meta.SCPCOA.TwistAng) * sind(meta.SCPCOA.GrazeAng))*180/pi;
    Multipath = mod(Azimuth - 180 + MultipathGround, 360);
    Shadow = mod(Azimuth - 180, 360);
end

%% Setup figure/axes
if isempty(p.Results.handle)
    f = figure;
    fig_pos = get(f,'Position');
    fig_pos(3) = 500;
    set(f,'Position',fig_pos);
    h = gca(f);   
else
    h = p.Results.handle;
end
axes(h);
delete(get(h,'Children'));
set(h,'Units','pixels');
axes_pos = get(h,'Position');
set(h,'Color',[0 0 0]);
set(h,'XLim',[0 190]);
set(h,'YLim',[0 210]);
set(h,'XTick',[]);
set(h,'YTick',[]);

% Text annotation
FontSize = floor(axes_pos(4)/20); %set font size based on height and width
text_flds = {'BackgroundColor',[0 0 0],'FontSize',FontSize};
white_flds = [text_flds {'Color','w'}];
text(8,193,strjoin(IID_line, ' / '),white_flds{:});
if exist('Tgt','var')
    geo_line = {['Tgt: ' strtrim(Tgt)]};
elseif exist('SCP','var')
    geo_line = {['Geo: ' latlonstr(Lat,'lat','include_symbols',false) ...
        '/' latlonstr(Lon,'lon','include_symbols',false)]};
else
    geo_line = {};
end
if exist('CC','var')
    geo_line = [geo_line ['CC: ' CC]];
end
text(8,171,strjoin(geo_line, ' / '),white_flds{:});
if exist('AzIPR','var') && exist('RgIPR','var')
    if ((AzIPR/RgIPR)-1) < .2
        res_line = {['IPR: ' num2str((AzIPR+RgIPR)/2, '%.1f') ' ft']};
    else %specify Az and Rg IPR separately
        res_line = {['IPR: ' num2str(AzIPR, '%.1f') '/' num2str(RgIPR, '%.1f') 'ft(A/R)']};
    end        
elseif exist('BW','var')
    res_line = {['BW: ' num2str(BW, '%.0f') ' MHz']};
else
    res_line = {};
end
if exist('RNIIRS','var')
    res_line = [res_line ['RNIIRS: ' RNIIRS]];
end
text(8,149,strjoin(res_line, ' / '),white_flds{:});
if exist('CollectDuration','var')
    cdp_line = {['CDP: ' num2str(CollectDuration, '%.1f') ' s']};
else
    cdp_line = {};
end
if ~isfield(meta,'CollectionInfo') || ~isfield(meta.CollectionInfo,'RadarMode') || ...
        ~isfield(meta.CollectionInfo.RadarMode,'ModeType') || ...
        ~strcmpi(meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    if exist('eff_ap_time','var')
        cdp_line = [cdp_line [num2str(eff_ap_time, '%.1f') ' s']];
    elseif isfield(meta, 'CollectionInfo') && ...
            isfield(meta.CollectionInfo, 'RadarMode') && ...
            isfield(meta.CollectionInfo.RadarMode,'ModeID')
        cdp_line = [cdp_line meta.CollectionInfo.RadarMode.ModeID];
    end
end
if exist('Pol','var')
    cdp_line = [cdp_line ['POL: ' Pol]];
end
text(8,127,strjoin(cdp_line, ' / '),white_flds{:});
if isfield(meta,'SCPCOA')
    text(8,105,['Azimuth: '  num2str(Azimuth, '%.1f') char(176)],white_flds{:});
    text(8,83,['Graze: ' num2str(meta.SCPCOA.GrazeAng, '%.1f') char(176)],white_flds{:});
    text(8,61,['Layover: ' num2str(Layover, '%.0f') char(176)],text_flds{:},'Color',[1 .66 0]);
    text(8,39,['Shadow: ' num2str(Shadow, '%.0f') char(176)],text_flds{:},'Color',[0 .65 1]);
    text(8,17,['Multipath: ' num2str(Multipath, '%.0f') char(176)],text_flds{:},'Color',[1 0 0]);
    text(100,17,meta.SCPCOA.SideOfTrack,text_flds{:},'Color','y');
    % Draw flight direction arrow
    if meta.SCPCOA.SideOfTrack=='R'
        arrow([120 17],[180 17],'Color',[1 1 0],'Width',2);
    else
        arrow([180 17],[120 17],'Color',[1 1 0],'Width',2);
    end

    if (~isfield(meta,'Grid') || strcmp(meta.Grid.ImagePlane,'SLANT')) && ~p.Results.GroundProject    
        Shadow = Azimuth - 180 - MultipathGround;
        Multipath = Azimuth - 180;
        Layover = Layover - MultipathGround;
    end

    % Arrow plot
    ArrowCenter = [145 72];
    ArrowLength = 45;
    Arrow_flds = {'Width', 1 ,'Length', 8, 'BaseAngle', 90, 'TipAngle', 20};

    Layover = 90-(Layover-Azimuth);
    arrow(ArrowCenter, [ArrowCenter(1)+ArrowLength*cosd(Layover) ArrowCenter(2)+ArrowLength*sind(Layover)], 'Color',[1 .66 0], Arrow_flds{:});

    % Shadow is always Down in Ground Plane
    Shadow = 90-(Shadow-Azimuth);
    arrow(ArrowCenter, [ArrowCenter(1)+ArrowLength*cosd(Shadow) ArrowCenter(2)+ArrowLength*sind(Shadow)], 'Color',[0 .65 1], Arrow_flds{:});

    % North Arrow
    North = Azimuth+90;
    arrow(ArrowCenter, [ArrowCenter(1)+(ArrowLength-5)*cosd(North) ArrowCenter(2)+(ArrowLength-5)*sind(North)], 'Color',[.58 .82 .31], Arrow_flds{:});
    text(175,33,'N','BackgroundColor',[0 0 0],'FontSize',FontSize-2,'Color',[.58 .82 .31]);

    % Multipath
    Multipath = North-Multipath;
    arrow(ArrowCenter, [ArrowCenter(1)+ArrowLength*cosd(Multipath) ArrowCenter(2)+ArrowLength*sind(Multipath)], 'Color',[1 0 0], Arrow_flds{:});
end

set(h,'Units','normalized');

end