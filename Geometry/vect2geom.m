function geometry = vect2geom( AIM, P, VEL, N )
%VECT2GEOM Compute SAR collection geometry from pointing vectors
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Conventions
%   - vectors are indicated by variables in caps
%   - unit vectors are denoted by a following underscore
%
% INPUTS:
%   AIM - required: Ground reference point in ECF
%   P   - required: Platform position, center aperture
%   VEL - required: Platform velocity
%   N   - required: Focus plane normal (defaults to WGS84 tangent plane
%                   unit normal)
%
% OUTPUTS:
%   geometry - SAR collection geometry structure
%
% VERSION:
%   1.0
%       -initial version
% TODO:

if nargin<4
    N = wgs_84_norm(AIM); % Defaults to WGS84 tangent plane unit normal
end

% range from aim point to platform
R  = P - AIM;
R_ = R / norm( R );
P_ = P / norm( P );

% aim point unit vectors
KDP_ = N / norm( N );
%KDP_ = AIM / norm( AIM );

JDP_ = proj( R_, KDP_ );
JDP_ = JDP_ / norm( JDP_ );

IDP_ = cross( JDP_, KDP_ );

% sensor azimuth
% azimuth = -atan2( IDP_(3), JDP_(3) );
ProjN_ = proj( [ 0; 0; 1 ], KDP_ );
ProjN_ = ProjN_ / norm( ProjN_ );

azimuth = atan2( dot( cross( JDP_, ProjN_ ), KDP_ ), dot( ProjN_, JDP_ ) );

if azimuth < 0
  azimuth = azimuth + 2 * pi;
end

%trajectory unit vector
TRAJ_ = VEL / norm( VEL );

% slant plane normal with ambiguous sense of up
SLANT = cross( R_, TRAJ_ );
SLANT_ = SLANT / norm( SLANT );

% direction of flight: sense < 0 - right, sense > 0 - left
sense = dot( SLANT_, KDP_ );
sense = sign( sense );

% corrected sense of up for slant plane normal
SLANT_ = sense * SLANT_;

% slope angle ( >= graze )
slope = acos( dot( SLANT_, KDP_ ) );

% grazing angle
graze = asin( dot( R_, KDP_ ) );

% squint angle
Vproj_ = proj( VEL, P );
Vproj_ = Vproj_ / norm( Vproj_ );
Rproj_ = proj( -R, P );
Rproj_ = Rproj_ / norm( Rproj_ );
squint = atan2( dot( cross( Vproj_, Rproj_ ), P_ ), dot( Rproj_, Vproj_ ) );

% layover angle
TMP_ = cross( SLANT_, KDP_ );
TMP_ = TMP_ / norm( TMP_ );
layover = asin( dot( JDP_, TMP_ ) );

% Doppler cone angle
dca = -sense * acos( dot( -R_, TRAJ_ ) );

% ground track angle in tangent plane at AIM point
TRACK_ = proj( TRAJ_, KDP_ );
TRACK_ = TRACK_ / norm( TRACK_ );
track = -sense * acos( dot( -1.*JDP_, TRACK_ ) );

% flight elevation
felev = asin( dot( KDP_, TRAJ_ ) );

% slant plane tilt angle
tilt = -acos( cos(slope) / cos(graze) ) * sign(layover);

% multipath angle
multipath = -atan( tan(tilt)*sin(graze) );

% results
geometry.azimuth   = azimuth;
geometry.graze     = graze;
geometry.slope     = slope;
geometry.squint    = squint;
geometry.layover   = layover;
geometry.multipath = multipath;
geometry.dca       = dca;
geometry.tilt      = tilt;
geometry.track     = track;
geometry.felev     = felev;

% left / right flight
if sense < 0
  geometry.right = 1;
elseif sense > 0
  geometry.right = -1;
else
  geometry.right = 0;
end

% ascending / descending
if VEL(3) > 0
  geometry.ascend = 1;
elseif VEL(3) < 0
  geometry.ascend = -1;
else
  geometry.ascend = 0;
end

end

% Project vector 'v' onto the plane defined by 'normal'
function p = proj( v, normal )
    n = normal / norm( normal );
    p = v - dot( v, n ) * n;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
