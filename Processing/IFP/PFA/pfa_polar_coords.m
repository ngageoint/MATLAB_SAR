function [k_a, k_sf] = pfa_polar_coords(pos, scp, ip_coa_pos, ipn, fpn)
%PFA_POLAR_COORDS Computes info required for mapping phase history to polar coordinates
%
% Inputs:
%     pos:       Nx3 array containing N positions of the sensor.
%     scp:       Scene center point.  Should be a 3-element vector.
%     coa_pos:   Position representing center of aperture or zero polar
%                angle.  This is the reference sensor position with respect
%                to which all the k_a (polar angles) are computed.  Should
%                be a 3-element vector.
%     ipn:       Image plane normal unit vector.  All positions are
%                projected into this plane before geometry information is
%                computed.  Should be a 3-element vector.
%     fpn:       Focus plane normal unit vector.  All projections are
%                performed orthoganally to this plane.  Should be a
%                3-element vector.
% Outputs:
%     k_a        Also known as "polar angle".  It is the angle (in radians)
%                in the image plane between each sensor position, the
%                center-of-aperture position, and the scene center point.
%                k_a is measured counter-clockwise.  k_a could roughly be
%                viewed as the "squint" angle for each pulse for a sensor
%                flying a straight-line.  (For a circular flight path,
%                squint would be constant, but this k_a angle will change.)
%     k_sf       The frequency scaling factor due to motion of the sensor
%                outside of the image plane.  Used to scale RF frequency to
%                radial spatial frequency in the image plane.
%
% Authors: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Handle input arguments
if all([numel(scp) numel(scp) numel(scp) numel(scp)]==3) && (size(pos,2)==3)
    % Assure correct vector orientations
    scp = scp(:)';
    ip_coa_pos = ip_coa_pos(:)';
    ipn = ipn(:)';
    fpn = fpn(:)';
else
    error('PFA_POLAR_COORDS:INVALID_INPUT','First argument POS should be Nx3.  The rest should be 3-element vectors.'); 
end

%% Project all positions to image plane along a line normal to the focus plane
ip_pos = proj(pos, fpn, scp, ipn);
ip_coa_pos = proj(ip_coa_pos, fpn, scp, ipn);

% Establish image plane coordinate system
ipx = ip_coa_pos - scp; % X direction is from scp toward coa (projected) position
ipx = ipx./sqrt(sum(ipx.^2)); % unit vector
ipy = cross(ipx,ipn); % Y direction forms a right handed coordinate system
% We already have z direction (ipz = ipn)

%% Compute polar angle of sensor position in image plane coordinate system
% Increasing counter clockwise-- which is totally weird, but that's the convention...
ip_range_vectors = bsxfun(@minus,ip_pos,scp);
k_a = -atan2(ip_range_vectors*ipy(:),ip_range_vectors*ipx(:));

%% Compute the spatial frequency scale factor due to projection to the image plane
% Collection geometry:
range_vectors = bsxfun(@minus,pos,scp);
range_vectors = bsxfun(@rdivide, range_vectors, sqrt(sum(range_vectors.^2,2))); % make unit vectors
sin_graze = range_vectors*fpn(:); % Dot product of fpn and all range unit vectors = cos of incidence angles = sin of graze angles
% Geometry of positions projected into the image plane:
ip_range_vectors = bsxfun(@rdivide, ip_range_vectors, sqrt(sum(ip_range_vectors.^2,2))); % make unit vectors
sin_graze_ip = ip_range_vectors*fpn(:); % cos of incidence angle = sin of graze angle
% k_sf = cos(collection_graze)./cos(image_plane_graze);
k_sf = sqrt(1-(sin_graze.^2))./sqrt(1-(sin_graze_ip.^2)); % cos = sqrt(1-sin^2)

end

%% Projection function
% Projection of a point along a given direction to a plane is just the
% intersection of the line defined by that point (l0) and direction (l) and
% the plane defined by a point in the plane (p0) and the normal (p):
% l0 + ((l0 - p0).p/(l.p))*l
% where . represents the dot product.
function pos_out = proj(points, line_direction, point_in_plane, plane_normal)
d = (bsxfun(@minus,point_in_plane, points) * plane_normal(:)) ./...
    (line_direction * plane_normal(:)); % Distance from point to plane in line_direction
pos_out = points + (d * line_direction);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////