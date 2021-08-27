function [phd_trimmed, cphd_nb, meta] = ExtractChipPhaseHistory(cphd_filename,SRP_new,img_support_size,channel)
% EXTRACTCHIPPHASEHISTORY extracts user-specified area from phase history
%   data = ExtractChipPhaseHistory(cphd_filename,img_offset,img_width,channel)
%   filters the input phase history to produce a reduced size phase history
%   that supports full resolution imaging of a user specified image chip
%   (slant plane rectangular area of dimensions IMG_SUPPORT_SIZE centered
%   on SRP_NEW.)
%
%   Separable 1-D range and 1-D azimuth processing, range first, azimuth
%   second, because of memory constraints.
%
% ASSUMPTIONS:
%   1) Assumes a spotlight collect, since it is assumed all points on the
%   ground are equally illuminated by all pulses.
%   2) Assumes a (nearly) monostatic collect.
%   3) This code assumes that the full phase history might not fit into
%   memory, but that the range trimmed phase history must be able to fit
%   into memory.
%   4) Also currently assumes that all pulses have common bandwidths and
%   frequency sampling.  If this is not the case, things become messy.
%   Azimuth filtering becomes difficult/impossible for pulses with
%   arbitrarily different frequency content.  Furthermore, its not clear
%   what the Fx0/F_SS of the resulting decimated pulses would be.
%
% INPUTS:
%   cphd_filename    : string : input phase history file
%   SRP_new          : [1x3]  : ECEF coordinates of new scene reference
%                               point.
%   img_support_size : If scalar, this is the radius of a sphere around the
%                      SRP to extract.
%                      If [1x2], this is the image chip extent in meters
%                      [azimuth_extent range_extent]  (in the slant plane).
%   channel          : scalar : Channel in the phase history file to chip.
%                               (Default = 1).
%
% OUTPUTS:
%   phd_trimmed : Phase history trimmed as requested with pulses in columns
%   cphd_nb     : CPHD formatted per pulse narrowband data for phd_trimmed
%   meta        : SICD-like meta structure
%
% Dan Hack, AFIT/ENG, 09AUG11
% Tim Cox, NRL, 19JUL12 - Modifications to bring to CPHD "X" format.
% Wade Schwartzkopf, NGA/IDT, 01DEC12 - Documentation and code clean up.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%%%%%%%%%%%%%%%%%%
% Constants      %
%%%%%%%%%%%%%%%%%%
c = 299792458; % Speed of light (m/s)
% Algorithm parameters
show_waitbar = true;
pBlockSize = 2000;
verification_pulse = 0; % Pulse number (or set to zero for no verification plots)

%%%%%%%%%%%%%%%%%%
% Pre-processing %
%%%%%%%%%%%%%%%%%%

% Default input parameters
if ~exist('channel','var')
    channel = 1; % Default
end
% Get the narrowband data  
reader_obj = open_ph_reader(cphd_filename);
meta = reader_obj.get_meta();
N = double(meta.Data.Channel(channel).NumVectors);
M = double(meta.Data.Channel(channel).NumSamples);  % num_samples per pulse
[ignore, all_nb] = reader_obj.read_cphd(1:N,[],1);
all_nb.ARP = (all_nb.TxPos + all_nb.RcvPos)/2; % Aperture reference position (nearly monostatic assumption)
% This function assumes spotlight data.
if isfield(meta, 'CollectionInfo') && isfield(meta.CollectionInfo, 'RadarMode') && ...
        isfield(meta.CollectionInfo.RadarMode, 'ModeType') && ...
        ~strcmpi(meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    error('EXTRACTCHIPPHASEHISTORY:NON_SPOTLIGHT_DATA',...
        'ExtractChipPhaseHistory function requires spotlight data.');
end
% Determine some collection parameters
[resolution, extent, delta_azimuth, total_azimuth] = pulse_info_to_resolution_extent(...
    all_nb.ARP([1 end],:) - all_nb.SRPPos([1 end],:), ... % Line-of-sight vector between ARP and SRP
    max(all_nb.SC0  + (all_nb.SCSS * (M-1))),... % Highest frequency sample (which has the least azimuth extent, so most conservative)
    max(all_nb.SCSS),... % Largest sample spacing (which has the least range extent)-- although we don't use it here
    [],... % This parameter should be bandwidth, but we don't need it here, since we aren't using resolution.
    N);

%%%%%%%%%%%%%%%%%%%%%%%
% 1:  Range Windowing %
%%%%%%%%%%%%%%%%%%%%%%%

% Slice out the required range bins from each pulse.

% First we determine how much range distance we need to cut out of each
% pulse.  
if isscalar(img_support_size) % A sphere around the SRP, the easiest.
    image_span_distance = img_support_size;
elseif numel(img_support_size)==2
    % Extracting range/azimuth oriented rectangle in the slant plane of
    % dimensions img_support_size ([azimuth, range]).
    %
    % In theory, since the geometry (and potentially frequency content) of
    % each pulse varies, each pulse could require a (potentially
    % drammatically) different set of range bins to exactly span our slant
    % plane rectangle. For the selected range bins to cover all corners of
    % the slant plane rectangle, including geometry and range curvature
    % contraints, the required distance span around the SRP for a single
    % pulse is:
    %    image_span_distance = 2*(sqrt(...
    %         (SRP_range*cos(theta) + (img_support_size(2)/2)).^2 + ...
    %         (SRP_range*abs(sin(theta)) + (img_support_size(1)/2)).^2) - ...
    %         SRP_range)
    % where
    %    theta is the angle of a given pulse off of the reference range
    %       direction for which our range/azimuth slant plane rectangle is
    %       defined.
    %    SRP_range is the distance from the ARP to SRP.
    %
    % Although the above equation may be exact per pulse, this distance
    % likely equates to a different number of range bins for each pulse,
    % and in MATLAB we prefer to work with rectangular arrays.  So for the
    % purposes of simplifying the computation, we pick the minimum distance
    % that spans our slant plane rectangle for all pulses.  This is not
    % optimum from the standpoint of the minimum amount of information
    % required to describe the requested image area, but it is required to
    % keep the computations tractable in MATLAB.
    % For simplicity, we only need to evaluate the above
    % image_span_distance equation with the pulse geometries at the ends of
    % the collect (when the theta angle, and thus the require range, is the
    % greatest) and with the minimum range to SRP for all pulses (since
    % that is when range curvature is at its greatest).
    min_SRP_range = min(sqrt(sum((all_nb.ARP - all_nb.SRPPos).^2,2)));
    image_span_distance = 2*(sqrt(...
        (min_SRP_range*cos(total_azimuth/2) + (img_support_size(2)/2)).^2 + ...
        (min_SRP_range*sin(total_azimuth/2) + (img_support_size(1)/2)).^2) - ...
        min_SRP_range);
    % Do we need some small amount of buffer on top of this?
else % For other areas (or volumes), this computation would have to change.
    error('EXTRACTCHIPPHASEHISTORY:UNRECOGNIZED_IMG_SUPPORT_SIZE',...
        'Unrecognized IMG_SUPPORT_SIZE format.');
end

% Determine range parameters
range_extent = c/(2*max(all_nb.SCSS));  % range profile extent (for pulse with the least extent)
rc_vec = fft_bin_pos_vec(M) * range_extent;  % Distance of each range bin from SRP
rc_ndx = find(abs(rc_vec) <= image_span_distance/2);  % range extraction indices
% These indices extract at least the range width needed to cover the
% requested area for all pulses.

% Allocate storage.  We assume this entire array of range-trimmed phase
% history can fit into memory.
M_new = length(rc_ndx); % Number of frequency samples after decimating
phd_trimmed_rangeonly = zeros(M_new,N);

% This is a vectorized implementation in which a block of pulses is
% processed on each iteration.  We process in blocks because it is likely
% that the entire untrimmed phase history cannot fit into memory.
%
% This loop processes all pulses in the phase history file because we
% assume a spotlight collect where all pulses equally affect our new
% reference point.  For non-spotlight collects, we should limit this loop
% to only pulses for which our new reference point fall in their beam.
% This would further reduce computation.

if show_waitbar
  wb = waitbar(0); tic
end

NumBlocks = ceil(N/pBlockSize); % Determine number of blocks.
for pbDx = 1:NumBlocks
    % Calculate pulse indices for this block
    pDx = ( ((pbDx-1)*pBlockSize+1):min(pbDx*pBlockSize,N) );
    
    % Read all samples from a block of pulses
    ph = reader_obj.read_cphd(pDx, 1:M, channel);
    
    % Calculate differential ranges (difference of range to current SRP and
    % range to desired SRP)
    delR = sqrt(sum((all_nb.ARP(pDx,:).' - repmat(SRP_new(:),[1 numel(pDx)])).^2)) - ...
           sqrt(sum((all_nb.ARP(pDx,:).' - all_nb.SRPPos(pDx,:).').^2));
        
    % Calculate frequency vector for each pulse: Fx0+Fx_SS*(0:M-1)
    f_orig = bsxfun(@plus,all_nb.SC0(pDx).',(0:M-1).'*all_nb.SCSS(pDx).');
    k_orig = 2*pi*f_orig/c;

    % Adjust phase reference (re-mocomp to new scene reference point)
    % This is currently just a linear shift based on difference in range
    % between the old and new reference points.  This is the 99% solution.
    % Might be some other higher-order components like tropospheric
    % adjustment.
    ph = ph.*exp(1i*2*bsxfun(@times,k_orig,delR));
    
    % Form range profile (IFFT)
    rc = ifft(ph,[],1);
    
    % Apply range gate and convert back to spatial freq domain
    phd_trimmed_rangeonly(:,pDx) = fft(rc(rc_ndx,:),[],1);

    % Generate plot comparisons for data from a single pulse
    if (verification_pulse ~= 0) && any(pDx == verification_pulse)
        % Recompute so we are using Fx_SS of this specific pulse, not max(Fx_SS)
        range_extent_orig = c/(2*all_nb.SCSS(verification_pulse));
        rc_vec_orig = fft_bin_pos_vec(M) * range_extent_orig;
        % Compute extent and range profile of trimmed pulse
        range_extent_new = c/(2*all_nb.SCSS(verification_pulse)*M/M_new);
        rc_vec_new = fft_bin_pos_vec(M_new) * range_extent_new;
        rc_new = ifft(phd_trimmed_rangeonly(:,verification_pulse));
        % Compare pre- and post- trimming
        figure;
        plot(fftshift(rc_vec_orig), fftshift(abs(rc(:,pDx==verification_pulse))), 'b', ...
            fftshift(rc_vec_new), fftshift(abs(rc_new)), 'r');
        set(gca,'xlim',[min(rc_vec) max(rc_vec)]);
        legend({'Original' 'Windowed'});
        xlabel('Range (m)'); ylabel('Magnitude');
        title(['Range Profile Comparison, Pulse ' num2str(verification_pulse)]);
    end
    
    if show_waitbar
        % Determine remaining execution time and display
        t_sofar = toc;
        t_est = (t_sofar*NumBlocks/pbDx)-t_sofar;
        wb_message=sprintf('Range Filtering, Time remaining: %s',datestr(datenum(0,0,0,0,0,t_est),13));
        waitbar(pbDx/NumBlocks,wb,wb_message);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2:  Azimuth (Doppler) Windowing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Due to re-mocomping to new reference point in the range processing step,
% the desired scene center is already at zero Doppler.
%
% Range filtering is exact.  Azimuth filtering (as implemented below) is
% approximate:
%
% 1) Pulses might not have uniform angular sampling. Since the FFT (and
% most filters) assume uniform sampling, this could be a point of concern.
% For now we assume that spacing is "close enough" to uniform, and that any
% distortions from the FFTing of non-uniform data will be somewhat undone
% with the inverse transform.  Also when we interpolate to new ARP
% positions later, we will retain this non-uniformity by interpolating with
% reference to pulse number rather than time.
% 2) Since we are filtering (at least roughly) across frequency sample
% here, there are some "polar formatting" approximations going on here.
%
% Note that if pulses have signficantly different Fx0 and/or Fx_SS, this
% step probably doesn't make any sense.

if show_waitbar
    waitbar(0,wb,'Azimuth Filtering');
end

% Azimuth profile vector
ac_vec = fft_bin_pos_vec(N)*extent(2); % extent(2) is azimuth extent

% Azimuth extraction indices
ac_ndx = find(abs(ac_vec) <= 1.1 * img_support_size(1)/2); % 10% extra buffer since our filtering is inexact 

% Number of azimuth frequency samples after decimating
N_new = length(ac_ndx);  

if show_waitbar
    waitbar(0.33,wb,'Azimuth Filtering');
end

% Azimuth compression (IFFT)
ac = ifft(phd_trimmed_rangeonly,[],2);

if show_waitbar
    waitbar(0.66,wb,'Azimuth Filtering');
end

% Window and azimuth uncompress (FFT)
phd_trimmed = fft(ac(:,ac_ndx),[],2);

if show_waitbar
    waitbar(1,wb,'Azimuth Filtering');
    pause(0.25);
    close(wb);
end

%%%%%%%%%%%%%%%%%%%%%%%
% Form pulse metadata %
%%%%%%%%%%%%%%%%%%%%%%%

% Generate output data structure.  Generic structure that reflects CPHD
% format metadata.

% Interpolate aperture reference position (ARP) -- assuming that
% interpolating the original trajectory in time will maintain phase
% coherency with respect to decimated pulses.
ARP_interp = interp1(1:N,all_nb.ARP,linspace(1,N,N_new),'spline');
% Time is really irrelevant in CPHD, but we compute it anyway.
time_ref = (all_nb.TxTime + all_nb.RcvTime)/2; % Monostatic approximation
time_interp = interp1(1:N,time_ref,linspace(1,N,N_new).','spline');

cphd_nb.TxTime = time_interp;
cphd_nb.TxPos = ARP_interp;
cphd_nb.RcvTime = time_interp;
cphd_nb.RcvPos = ARP_interp;
cphd_nb.SRPPos = repmat(SRP_new(:).',[N_new 1]);
% Interpolation along Fx probably isn't valid if Fx values jitter.
cphd_nb.SC0 = interp1(1:N,all_nb.SC0,linspace(1,N,N_new).','spline');
cphd_nb.SCSS = interp1(1:N,all_nb.SCSS*(M/M_new),linspace(1,N,N_new).','spline');
cphd_nb.FX1 = interp1(1:N,all_nb.FX1,linspace(1,N,N_new).','spline');
cphd_nb.FX2 = interp1(1:N,all_nb.FX2,linspace(1,N,N_new).','spline');

end

% Compute a vector of the positions of the bins in an FFT'd vector.
% Output is in normalized units (1 is Nyquist).  Need to multiply by extent
% to convert to real units like meters (or Hz).
function bin_pos_norm = fft_bin_pos_vec(Nfft)
    bin_pos_norm = (0:(Nfft-1))/Nfft; % Normalized units (no sampling frequency used)
    ind_to_flip = (floor((Nfft-1)/2)+2):Nfft; % These positions are 0.5 or above.  Should be negative units.
    bin_pos_norm(ind_to_flip) = bin_pos_norm(ind_to_flip) - 1;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////