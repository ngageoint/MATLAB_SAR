function [ sicd_meta ] = meta2sicd_palsar2ledimg( led_meta, img_meta )
%META2SICD_PALSAR2LEDIMG Computes SICD parameters that need information
% from both the IMG and LED files in the ALOS PALSAR 2 package
%
% Takes as input metadata structures from read_ceos_led_meta and read_ceos_img_meta.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% ImageData
% We choose the natively defined ctr_pixel/ctr_line as the SCP, but this
% could really be set to anything.
sicd_meta.ImageData.SCPPixel.Row = uint32(led_meta.data.ctr_pixel);
% ctr_line in native metadata is measured from early edge of image, whereas
% SICD SCP is measured from left of image when in a view-from-above
% orientation.
if led_meta.data.clock_angle<0  % Left looking
    sicd_meta.ImageData.SCPPixel.Col = uint32(...
        (img_meta.num_lines-1) - led_meta.data.ctr_line);
else
    sicd_meta.ImageData.SCPPixel.Col = uint32(led_meta.data.ctr_line);
end

%% Position
% We make the time of the first line as Timeline.CollectStart, and thus the
% reference time for the SICD, so position polynomials must be in relation
% to it.
ref_time_offset = img_meta.signal(1).usec/1e6 - led_meta.pos.sec;
if led_meta.pos.day_in_year ~= img_meta.signal(1).day  % Assume we crossed midnight
    ref_time_offset = ref_time_offset + (24 * 60 * 60);  % So add a day of seconds
end
% Get state vectors
state_vector_T = ((0:(led_meta.pos.num_pts-1)) * led_meta.pos.int).' - ...
    ref_time_offset;
state_vector_pos = zeros(led_meta.pos.num_pts,3);
state_vector_vel = zeros(led_meta.pos.num_pts,3);
for i=1:led_meta.pos.num_pts
    state_vector_pos(i,1) = led_meta.pos.pts(i).pos_x;
    state_vector_pos(i,2) = led_meta.pos.pts(i).pos_y;
    state_vector_pos(i,3) = led_meta.pos.pts(i).pos_z;
    state_vector_vel(i,1) = led_meta.pos.pts(i).vel_x;
    state_vector_vel(i,2) = led_meta.pos.pts(i).vel_y;
    state_vector_vel(i,3) = led_meta.pos.pts(i).vel_z;
end
% Limit state vectors so polynomial fit is only performed for times
% moderatley close to this collect.  PALSAR state vectors are typically
% spaced fairly far apart (usually 60 seconds), so we don't get that many
% actually inside the collect (often only 1).  We use +/- 4 minutes around
% collect, so this should at least include 9 state vectors, which should
% support fitting a 6th order polynomial and result in precisions of about
% a mm.
buffer = 60*4;
duration = (img_meta.signal(end).usec - img_meta.signal(1).usec) * 1e-6;
valid = (state_vector_T + buffer > 0) & ...
    (state_vector_T - duration - buffer < 0);  % Within collect +/- buffer
% sv2poly.m shows ways to determine best polynomial order.  Usually 5th
% order is sufficient to accurately describe any typical spaceborne SAR
% collect. However, in this case, since the polynomial fit is over such a
% large period of time, we bump up the polynomial order to 6.
polyorder = min(6, numel(state_vector_T) - 1);
old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
P_x = polyfit(state_vector_T(valid), state_vector_pos(valid,1), polyorder);
P_y = polyfit(state_vector_T(valid), state_vector_pos(valid,2), polyorder);
P_z = polyfit(state_vector_T(valid), state_vector_pos(valid,3), polyorder);
warning(old_state);
sicd_meta.Position.ARPPoly.X = P_x(end:-1:1).';
sicd_meta.Position.ARPPoly.Y = P_y(end:-1:1).';
sicd_meta.Position.ARPPoly.Z = P_z(end:-1:1).';

%% RMA
% These Grid fields are redundant from meta2sicd_palsar2led, but are needed
% extensively in the RMA section.
sicd_meta.Grid.Row.SS = led_meta.data.pixel_spacing;
sicd_meta.Grid.Col.SS = led_meta.data.line_spacing;
% Zero doppler spacing in seconds.  Could also use line timing in img file.
ss_zd_s = 1000/led_meta.data.prf;
if led_meta.data.clock_angle<0  % Left looking
    % scp_az_offset is the distance (in pixels) from the SCP to the left
    % hand side of the image when in a "view-from-above" orientation, as is
    % used is SICD.
    scp_az_offset = double((img_meta.num_lines-1) - sicd_meta.ImageData.SCPPixel.Col);
    ss_zd_s = -ss_zd_s;
else
    scp_az_offset = double(sicd_meta.ImageData.SCPPixel.Col);
end
sicd_meta.RMA.INCA.TimeCAPoly = [abs(ss_zd_s) * scp_az_offset ; ... % SCP
    ss_zd_s/led_meta.data.line_spacing];  % Convert spacing from sec/pixels to sec/meters
sicd_meta.RMA.INCA.R_CA_SCP = img_meta.signal(1).slant_rng + ...
    sicd_meta.Grid.Row.SS * double(sicd_meta.ImageData.SCPPixel.Row);
% DopCentroidPoly
dop_poly_az=[led_meta.data.at_dop_const, led_meta.data.at_dop_lin, led_meta.data.at_dop_quad];
dop_poly_rg=[led_meta.data.xt_dop_const, led_meta.data.xt_dop_lin, led_meta.data.xt_dop_quad];
sicd_meta.RMA.INCA.DopCentroidPoly=zeros(3);
sicd_meta.RMA.INCA.DopCentroidPoly(1) = ...  % Compute doppler centroid value at SCP
    polyval(dop_poly_rg(end:-1:1), double(sicd_meta.ImageData.SCPPixel.Row)) + ...
    polyval(dop_poly_az(end:-1:1), scp_az_offset) - ...
    mean([dop_poly_az(1) dop_poly_rg(1)]); % These should be identical
% Convert polynomials from Hz/pixel from edge to Hz/meter from SCP
dop_poly_az_shifted=polyshift(dop_poly_az, scp_az_offset);
dop_poly_rg_shifted=polyshift(dop_poly_rg, double(sicd_meta.ImageData.SCPPixel.Row));
% Scale 1D polynomials to from Hz/pixel to Hz/m^n
dop_poly_az_scaled=dop_poly_az_shifted.*...
    sign(ss_zd_s).*sicd_meta.Grid.Col.SS.^(0:(numel(dop_poly_az)-1));
dop_poly_rg_scaled=dop_poly_rg_shifted.'.*...
    sicd_meta.Grid.Row.SS.^(0:(numel(dop_poly_rg)-1)).';
sicd_meta.RMA.INCA.DopCentroidPoly(2:end,1)=dop_poly_rg_scaled(2:end);
sicd_meta.RMA.INCA.DopCentroidPoly(1,2:end)=dop_poly_az_scaled(2:end);
sicd_meta.RMA.INCA.DopCentroidCOA=true;
sicd_meta.Grid.Col.DeltaKCOAPoly=...
    sicd_meta.RMA.INCA.DopCentroidPoly*ss_zd_s/sicd_meta.Grid.Col.SS;
% Compute DRateSFPoly
dop_rate_poly_rg=[led_meta.data.xt_dop_rate_const, ...
    led_meta.data.xt_dop_rate_lin, led_meta.data.xt_dop_rate_quad];
dop_rate_poly_rg_shifted=polyshift(dop_rate_poly_rg, ...
    double(sicd_meta.ImageData.SCPPixel.Row));
dop_rate_poly_rg_scaled=dop_rate_poly_rg_shifted.'.*...
    sicd_meta.Grid.Row.SS.^(0:(length(dop_rate_poly_rg)-1)).';
% For the purposes of the DRateSFPoly computation, we ignore any
% changes in velocity or doppler rate over the azimuth dimension.  These
% are generally insignificant anyway.
vel_x = polyval(polyder(P_x), sicd_meta.RMA.INCA.TimeCAPoly(1));
vel_y = polyval(polyder(P_y), sicd_meta.RMA.INCA.TimeCAPoly(1));
vel_z = polyval(polyder(P_z), sicd_meta.RMA.INCA.TimeCAPoly(1));
vm_ca_sq = vel_x.^2 + vel_y.^2 + vel_z.^2; % Magnitude of the velocity squared
r_ca = [sicd_meta.RMA.INCA.R_CA_SCP; 1]; % Polynomial representing range as a function of range distance from SCP
fc = SPEED_OF_LIGHT() / led_meta.data.wavelength;
sicd_meta.RMA.INCA.DRateSFPoly = - conv(dop_rate_poly_rg_scaled,r_ca) * ... % Multiplication of two polynomials is just a convolution of their coefficients
    SPEED_OF_LIGHT / (2 * fc * vm_ca_sq(1)); % Assumes a SGN of -1
% TimeCOAPoly
% TimeCOAPoly=TimeCA+(DopCentroid/dop_rate)
% Since we can't evaluate this equation analytically, we will evaluate
% samples of it across our image and fit a 2D polynomial to it.
POLY_ORDER = 2; % Order of polynomial which we want to compute
grid_samples = POLY_ORDER + 1; % in each dimension
coords_az_m = linspace(-double(sicd_meta.ImageData.SCPPixel.Col),...
    double(img_meta.num_lines-sicd_meta.ImageData.SCPPixel.Col-1), grid_samples) * ...
    sicd_meta.Grid.Col.SS;
coords_rg_m = linspace(-double(sicd_meta.ImageData.SCPPixel.Row),...
    double(img_meta.num_pixels-sicd_meta.ImageData.SCPPixel.Row-1), grid_samples) * ...
    sicd_meta.Grid.Row.SS;
timeca_sampled = sicd_polyval2d(sicd_meta.RMA.INCA.TimeCAPoly(:).',coords_az_m,coords_rg_m);
dopcentroid_sampled = sicd_polyval2d(sicd_meta.RMA.INCA.DopCentroidPoly,coords_az_m,coords_rg_m);
doprate_sampled = sicd_polyval2d(dop_rate_poly_rg_scaled.',coords_az_m,coords_rg_m);
timecoapoly_sampled = timeca_sampled+(dopcentroid_sampled./doprate_sampled);
% Least squares fit for 2D polynomial
% A*x = b
[coords_az_m, coords_rg_m] = ndgrid(coords_az_m, coords_rg_m);
a = zeros(grid_samples^2, (POLY_ORDER+1)^2);
for k = 0:POLY_ORDER
    for j = 0:POLY_ORDER
        a(:,k*(POLY_ORDER+1)+j+1) = (coords_rg_m(:).^j).*(coords_az_m(:).^k);
    end
end
b_coa = zeros((POLY_ORDER+1)^2,1);
for k=1:((POLY_ORDER+1)^2)
   b_coa(k)=sum(timecoapoly_sampled(:).*a(:,k)); % center of aperture
end
A=zeros((POLY_ORDER+1)^2);
for k=1:((POLY_ORDER+1)^2)
    for j=1:((POLY_ORDER+1)^2)
        A(k,j)=sum(a(:,k).*a(:,j));
    end
end
old_warning_state=warning('off','MATLAB:nearlySingularMatrix');
x=A\b_coa; % MATLAB often flags this as badly scaled, but results still appear valid
warning(old_warning_state);
sicd_meta.Grid.TimeCOAPoly=reshape(x, POLY_ORDER+1, POLY_ORDER+1);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////