function Az = ComputeAz(meta, percent, pixel_coords)
% COMPUTAZ Compute azimuth angle (angle from north) across spatial frequency
%
% Should work generically for many types of SAR complex data, not just
% spotlight collects.
%
% INPUTS:
%   meta         - required: SICD metadata structure
%   percent      - optional: Array of values from 0 to 100.  Indicates the
%                  fraction across the spatial frequency in azimuth for
%                  which to compute azimuth angle (not including zeropad).
%                  Default is scene center point center of aperture.
%   pixel_coords - optional: [column_index, row_index] (az,rng) index to
%                  point of interest.  Default is scene center point (SCP).
%
% OUTPUTS:
%   Az - computed azimuth angle for each PERCENT value
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Handle default values for input parameters
% Default for percent is just SCPCOA azimuth
% If percent is not passed, we might not need TimeCOAPoly, so we check this first.
if ~exist('percent','var')
    % SCPCOA.AzimAng is an actual field in latest SICD spec
    if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'AzimAng')
        Az = meta.SCPCOA.AzimAng;
        return;
    end
    percent = 50; % If SCPCOA.AzimAng not available, we will compute it.
end
if ~isfield(meta.Grid,'TimeCOAPoly') % TimeCOAPoly field is required
    % For spotlight, we can take some shortcuts with missing metadata
    if isfield(meta,'CollectionInfo') &&...
            isfield(meta.CollectionInfo,'RadarMode')&&...
            isfield(meta.CollectionInfo.RadarMode,'ModeType')&&...
            strcmpi(meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
        if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SCPTime') &&...
                isfield(meta,'Position') && isfield(meta.Position,'ARPPoly')
            meta.Grid.TimeCOAPoly = meta.SCPCOA.SCPTime;
        elseif isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'ARPPos') &&...
                isfield(meta.SCPCOA,'ARPVel')
            meta.Grid.TimeCOAPoly = 0;
            meta.Position.ARPPoly.X = [meta.SCPCOA.ARPPos.X; meta.SCPCOA.ARPVel.X];
            meta.Position.ARPPoly.Y = [meta.SCPCOA.ARPPos.Y; meta.SCPCOA.ARPVel.Y];
            meta.Position.ARPPoly.Z = [meta.SCPCOA.ARPPos.Z; meta.SCPCOA.ARPVel.Z];
        end
    end
    if ~isfield(meta.Grid,'TimeCOAPoly')
        error('COMPUTEAZ:TIMECOAPOLY_MISSING','Unable to compute SICD Grid.TimeCOAPoly field from complex data.');
    end
end
% Point of interest (POI) only required for spatially variant COA
if exist('pixel_coords','var') && ~isempty(pixel_coords)
    % Convention for point_slant_to_ground() pixel coordinates is reverse
    % of what is passed to this function (column/row), so we have to swap
    % values.    
    poi_ecf = geodetic_to_ecf(point_slant_to_ground([pixel_coords(2); pixel_coords(1)], meta)).';
    % Calculate center of aperture (COA) for given point
    time_coa = sicd_polyval2d(meta.Grid.TimeCOAPoly,pixel_coords(1),pixel_coords(2),meta);
else % Default to SCP
    poi_ecf = [meta.GeoData.SCP.ECF.X, meta.GeoData.SCP.ECF.Y, meta.GeoData.SCP.ECF.Z];
    time_coa = meta.Grid.TimeCOAPoly(1); % SCP COA time is this by definition
    if any(meta.Grid.TimeCOAPoly(2:end)~=0) % Spatially variant COA
        warning('COMPUTEAZ:PIXEL_COORDS_REQUIRED','For non-spotlight data, point of interest must be specified for accurate azimuth angles.');
    end
end

%% Calculate geometry info for center of aperture
pos_coefs=[meta.Position.ARPPoly.X(end:-1:1)...
    meta.Position.ARPPoly.Y(end:-1:1)...
    meta.Position.ARPPoly.Z(end:-1:1)]; % Position polynomial
ARP=[polyval(pos_coefs(:,1),time_coa)...
    polyval(pos_coefs(:,2),time_coa)...
    polyval(pos_coefs(:,3),time_coa)]; % Aperture position at COA
% Velocity polynomial is derivative of position polynomial
vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
ARV=[polyval(vel_coefs(:,1), time_coa)...
    polyval(vel_coefs(:,2), time_coa)...
    polyval(vel_coefs(:,3), time_coa)]; % Aperture velocity at COA
LOS = poi_ecf - ARP; % Line of site vector at COA
R = norm(LOS); % Range from sensor to POI at COA
V = norm(ARV); % Speed at COA
DCA =  acos((ARV/V)*(LOS.'/R)); % Doppler Cone Angle to POI at COA

%% Compute effective aperture positions for each percentage
Theta = meta.Grid.Col.ImpRespBW/meta.Grid.Row.KCtr; % Angular aperture spanned by this complex image
EffectiveDuration = (Theta*R/V)/sin(DCA); % Theta is expressed as a ratio, rather than an angle
% For the purposes of this computation, we assume a linear flight path
% centered around the COA ARP and ARV, and that the center of the deskewed
% spatial frequency in azimuth is the azimuth angle at COA.
t = EffectiveDuration * ((-1/2) + (percent(:) / 100));
pos = [ARP(1) + (t * ARV(1)), ...
    ARP(2) + (t * ARV(2)), ...
    ARP(3) + (t * ARV(3))];
% We could also compute position given the position given the full position
% polynomial.  However, given that the PERCENT parameter is actually a
% fraction of the spatial frequency space and not time, it is not clear
% that this will be any more accurate.
% t = time_coa + (EffectiveDuration * ((-1/2) + (percent(:) / 100)));
% pos=[polyval(pos_coefs(:,1),t)...
%     polyval(pos_coefs(:,2),t)...
%     polyval(pos_coefs(:,3),t)];

%% Calculate actual azimuth angles
gpn = wgs_84_norm(poi_ecf); % Ground plane normal, based on wgs_84 ellipsoid
range = pos-repmat(poi_ecf,[size(t,1),1]); % Range vectors for each time
north_ground = [ 0 0 1 ]-([ 0 0 1 ]*gpn)*gpn.'; % Project north onto ground plane
range_ground=range-(range*gpn)*gpn.'; % Project range vector onto ground plane
Az = atan2( cross( range_ground, repmat(north_ground,[size(t,1),1]), 2 ) * gpn, ...
    range_ground * north_ground.')*180/pi; % Angle between the two projected vectors
Az(Az < 0) = Az(Az < 0) + 360; % Assure in [0,360], not [-180,180]
Az = reshape(Az,size(percent)); % Reform to shape of original input

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////