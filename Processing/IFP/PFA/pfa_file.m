function [sicd_meta] = pfa_file(input_filename_cphd, output_filename_sicd, varargin)
%PFA_FILE Implements the polar format IFP algorithm on a file of
%motion-compensated phase history data in CPHD (version 3.0) file format
%
% This function actually contains almost no actual PFA code directly, but
% rather calls PFA helper functions that do all of the real PFA work.  This
% function does however demonstrates how to break up a file and process it
% in sections through intermediate files.
%
% PFA_FILE(INPUT_FILENAME_CPHD, OUTPUT_FILENAME_SICD, 'PropertyName', PropertyValue, ...)
%
%       Property name     Description
%       resolution        Desired resolution 3dB IPR width
%                         ([range_resolution azimuth_resolution]) in
%                         meters.  The necessary pulses and samples for
%                         image formation at this resolution will be
%                         calculated (taken from middle of data.) If
%                         pulse_range and/or sample_range properties are
%                         defined, then this parameter is ignored.  Default
%                         is the full resolution that the collect supports.
%       pulse_range       Set of pulses to use for image formation. (For
%                         example, 1:1000.) Default is to calculate this
%                         from resolution property.
%       sample_range      Set of samples to use for image formation. (For
%                         example, 1:1000.) Default is to calculate this
%                         from resolution property.
%       channel           Channel to use from CPHD file.  Default is 1.
%       sample_rate       Samples per IPR.  Default is 1.5.
%       max_block_size    When data can't be processed within memory, this
%                         is the size of the blocks (in bytes) to use
%                         process the data.
%       quiet             If false, this reports stats on collection and
%                         IFP parameters.  Default is true.
%
% Note for larger datasets, this routine executes "out of core" storing of
% intermediate results in temporary files, which should be automatically
% cleaned up upon completion of this function.  However, if function is
% stopped in the middle, either through an error or user intervention, some
% of the files may remain in MATLAB's temporary directory (tempdir).
% 
% TODO:
%   - Check bistatic frequency scaling
%   - More input parameters (ipn, fpn, coa, ss, etc.)
%   - Use sinc interpolation instead of linear interpolation
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters and read metadata
ph_reader = open_ph_reader(input_filename_cphd);
ifp_params = select_pulses_samples_cphd(ph_reader, varargin{:}); % Parse generic IFP parameters
p1 = inputParser; % Parse PFA-specific parameters
p1.KeepUnmatched=true;
p1.addParamValue('max_block_size',[], @isscalar);
p1.addParamValue('sample_rate',1.5, @isscalar);
p1.FunctionName = mfilename;
p1.parse(varargin{:});
ifp_params.sample_rate = p1.Results.sample_rate;
max_block_size = p1.Results.max_block_size;
if isempty(max_block_size)
    if ispc
        [uv, sv] = memory;
        % max_block_size is just the size of the largest single array we
        % will handle at once.  We will need significantly more memory than
        % this.  A factor of 16 seems to be reasonable.  Its a nice round
        % power of two that seems to avoid swapping in limited testing.
        max_block_size = min(sv.PhysicalMemory.Available/16, ...
            uv.MaxPossibleArrayBytes);
    else % Can't gauge memory in UNIX, so just put arbitrary value
        max_block_size = 2^27; % About 100Meg
    end
end
cphd_meta = ph_reader.get_meta();
[ignore, nbdata] = ph_reader.read_cphd(ifp_params.pulse_range, [], 1);
if ~strcmpi(cphd_meta.CollectionID.RadarMode.ModeType,'SPOTLIGHT') || ...
        any(any(diff(nbdata.SRPPos)))  % Assure spotlight data
    error('PFA_FILE:UNSUPPORTED_COLLECT_TYPE','Unsupported collection mode.  Currently only spotlight data is supported.');
else
    scp = nbdata.SRPPos(1,:);
end
if ~strcmp(cphd_meta.Global.DomainType, 'FX')
    error('PFA_FILE:UNSUPPORTED_FILE_TYPE','Unsupported CPHD type.  Currently only FX domain is supported.');
end

%% Compute some values we will need
% Number of pulses and samples that we will be processing in this function.
% Defined here purely for conciseness of notation.
num_pulses  = length(ifp_params.pulse_range);
num_samples = length(ifp_params.sample_range);

% Decide whether we can process in memory, or whether we need to use
% intermediate files.  We hope to avoid using intermediate files since file
% reads and writes take the great majority of the processing time here.  In
% one profiling test, file reads and writes took 98% of the time of this
% funcion.  For this reason, we have attempted to minimize the number of
% reads and writes, and we try to comment in the code sections that affect
% these file I/O delays.  Unfortunately virtual memory and swapping (which
% is non-optimal file I/O out of our control) is much worse than the file
% I/O explicitly handled in this function, so try to set max_block_size to
% something as large as possible without the operating system resorting to
% virtual memory.
data_element_size = 8; % In bytes.  Assumes single precision data
process_in_blocks = max_block_size<(num_pulses*num_samples*prod(ifp_params.sample_rate)*data_element_size);
if process_in_blocks
    warning('PFA_FILE:FILE_PROCESSING','This data must be processed through intermediate files, rather than all in memory.  This could take a while...');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We'll generate the pulse sample points in a U,V coordinate system in
% which the V axis is parallel to the center of aperture and directed
% 'outward' (increasing frequency) in the k_x, k_y plane.  U is orthoganal
% to V and points left when viewed from above.
%
% In this function, (k_r, k_a) are the polar coordinates (radial and
% angular) of every pulse/sample, and (k_u, k_v) are the rectangular
% coordinates in the U,V coordinate system.  We will make the units of our
% k-space to be cycles/meter, which is consistent with SICD metadata.

[ bi_pos, bi_freq_scale ] = pfa_bistatic_pos( nbdata.TxPos, nbdata.RcvPos, nbdata.SRPPos );
% Temporary for testing fpn/ipn.  Actually this should be an input argument.
ifp_params.fpn=wgs_84_norm(scp).'; % Compute normal to WGS_84 ellipsoid
scp = [0 0 0]; % bistatic position is with respect to scp as origin
ifp_params.ref_pulse_index  = floor(num_pulses/2); % Center pulse seems as good as any...
arp_coa_vel = diff(bi_pos(ifp_params.ref_pulse_index + [1 -1],:));
arp_coa = bi_pos(ifp_params.ref_pulse_index,:);
srv=(arp_coa-scp).'; % slant range vector
look=ifp_params.fpn*cross(srv,arp_coa_vel).';
ifp_params.ipn=look*cross(srv,arp_coa_vel);
ifp_params.ipn=ifp_params.ipn/norm(ifp_params.ipn); % Slant plane unit normal
% end temporary portion here

[k_a, k_sf] = pfa_polar_coords(bi_pos, scp, arp_coa, ifp_params.ipn, ifp_params.fpn); % Angular coordinate of each pulse
% Converting from raw RF frequency of the received pulse to radial position
% of each pulse/sample in the image formation plane in "K-space" involves a
% few scaling factors:
rf_to_rad = (2/SPEED_OF_LIGHT) .* ... % Convert from cycles/second to cycles/meter
    k_sf .* ... % Compensate for out-of-plane motion by projecting into image formation plane
    bi_freq_scale; % For bistatic collects, a factor to account for using the equivalent monostatic position
k_r0 = rf_to_rad .* (nbdata.SC0 + (nbdata.SCSS*min(double(ifp_params.sample_range)-1))); % Radial position of the first sample in each pulse
k_r_ss = rf_to_rad .* nbdata.SCSS; % Radial position spacings between the samples in each pulse
% Compute new coordinates onto which to interpolate
[k_v_bounds, k_u_bounds] = pfa_inscribed_rectangle_coords(k_a, k_r0, ...
    double(max(ifp_params.sample_range) - min(ifp_params.sample_range)) * ...
    k_r_ss); % Effective bandwidth in image plane

% Compute the sicd metadata parameters
ifp_params.k_u_bounds = k_u_bounds;
ifp_params.k_v_bounds = k_v_bounds;
ifp_params.image_size_pixels = floor(ifp_params.sample_rate .* [num_samples num_pulses]);
ifp_params.k_a = k_a;
ifp_params.k_sf = k_sf;
sicd_meta = pfa_sicd_meta(cphd_meta,nbdata,ifp_params);
% Coordinates onto which we will interpolate
new_v = linspace(k_v_bounds(1), k_v_bounds(2), num_samples).';
new_u = linspace(k_u_bounds(1), k_u_bounds(2), num_pulses).';

%% Range interpolation across samples
% The interpolation will be done in two parts.  First we'll interpolate
% radially within each pulse to a constant spacing in V, then we'll use the
% just-interpolated values to interpolate across pulses (at each sample) to
% a constant spacing in U.

% Setup temporary file
% The range interpolation stage includes an in-order read (CPHD is written
% in "pulse-major" order) and an in-order write (data interpolated in
% range).  We have to "corner-turn" (or transpose) after range
% interpolation (as we do on all 4 steps), but it is typically faster to
% read from a file out-of-order than write out-of-order, so we will do the
% corner-turn in the read of the U-interpolation stage, rather than the
% write of this range interpolation stage.
if process_in_blocks
    filename_range_interp = tempname;
    tempsicdmeta.ImageData.NumRows   = num_pulses;
    tempsicdmeta.ImageData.NumCols   = num_samples;
    writer_range_interp = SIOWriter(filename_range_interp, tempsicdmeta, 'include_sicd_metadata', false);
end

% Iterate through sets of pulses for range interpolation
first_pulse_in_set=1; % Out of all pulses selected for processing
max_num_lines = floor(max_block_size/(num_samples*data_element_size));
wb_hand=waitbar(0,'V interpolation...');
while(first_pulse_in_set<num_pulses)
    waitbar((first_pulse_in_set-1)/num_pulses,wb_hand);
    pulse_indices = first_pulse_in_set:...
        min(num_pulses,first_pulse_in_set+max_num_lines-1);
    data_block = ph_reader.read_cphd(ifp_params.pulse_range(pulse_indices),...
        ifp_params.sample_range, ifp_params.channel);

    data_block = pfa_interp_range(data_block, k_a(pulse_indices), ...
        k_r0(pulse_indices),k_r_ss(pulse_indices), new_v);

    if process_in_blocks
        writer_range_interp.write_chip(data_block,[1 first_pulse_in_set]); % Save range interpolated data in this block
    end
    first_pulse_in_set=pulse_indices(end)+1;
end
if process_in_blocks
    clear writer_range_interp % Close file
end
ph_reader.close();

%% Interpolate U values across pulses
% Setup temporary files
if process_in_blocks
    reader_range_interp = open_reader(filename_range_interp);
    filename_u_interp = tempname;
    tempsicdmeta.ImageData.NumRows   = num_samples;
    tempsicdmeta.ImageData.NumCols   = num_pulses;
    writer_u_interp = SIOWriter(filename_u_interp, tempsicdmeta, 'include_sicd_metadata', false);
end

% Iterate through blocks of data for U interpolation
first_sample_in_block=1; % Out of all samples selected for processing
max_num_lines = floor(max_block_size/(num_pulses*data_element_size));
waitbar(0,wb_hand,'U interpolation...');
while(first_sample_in_block<num_samples)
    waitbar((first_sample_in_block-1)/num_samples,wb_hand);
    sample_range = [first_sample_in_block ...
        min(num_samples,first_sample_in_block+max_num_lines-1)];
    if process_in_blocks
        data_block = reader_range_interp.read_chip(sample_range, [1 num_pulses]);
    end
    % We transpose this data so we can process it along the first
    % dimensions, which is somewhat faster, and write it out in order,
    % which can be much faster.
    data_block = pfa_interp_azimuth(data_block.', k_a, ...
        new_v(sample_range(1):sample_range(2)), new_u);
    
    if process_in_blocks
        writer_u_interp.write_chip(data_block,[1 first_sample_in_block]);
    end
    first_sample_in_block=sample_range(2)+1;
end
if process_in_blocks
    clear writer_u_interp; % Close file
    reader_range_interp.close(); % Close file
    delete(filename_range_interp); % Intermediate file no longer needed
end

%% FFT V
% Setup temporary files
if process_in_blocks
    reader_u_interp = open_reader(filename_u_interp);
    filename_v_fft = tempname;
    tempsicdmeta.ImageData.NumRows   = num_pulses;
    tempsicdmeta.ImageData.NumCols   = ifp_params.image_size_pixels(1);
    writer_u_fft = SIOWriter(filename_v_fft, tempsicdmeta, 'include_sicd_metadata', false);
end

% Iterate through blocks of data one more time for V FFT
first_line_in_block=1; % Out of all lines along V in the resampled image
max_num_lines = floor(max_block_size/(length(new_v)*data_element_size));
waitbar(0,wb_hand,'V FFT...');
while(first_line_in_block<length(new_u))
    waitbar((first_line_in_block-1)/single(length(new_u)),wb_hand);
    line_range = [first_line_in_block ...
        min(length(new_u),first_line_in_block+max_num_lines-1)];
    if process_in_blocks
        data_block = reader_u_interp.read_chip(line_range, [1 length(new_v)]);
    end

    data_block = pfa_fft_zeropad_1d(data_block.', ifp_params.sample_rate);
    
    if process_in_blocks
        writer_u_fft.write_chip(data_block,[1 first_line_in_block]);
    end
    first_line_in_block=line_range(2)+1;
end
if process_in_blocks
    clear writer_u_fft % Close file
    reader_u_interp.close(); % Close file
    delete(filename_u_interp); % Intermediate file no longer needed
end

%% FFT U
if process_in_blocks
    reader_v_fft = open_reader(filename_v_fft);
end
writer_final_sicd = SICDWriter(output_filename_sicd, sicd_meta); % Final image will be written here

% Iterate through blocks of data for U FFT
first_line_in_block=1; % Out of all rows in the final image
max_num_lines = floor(max_block_size/(length(new_u)*data_element_size));
waitbar(0,wb_hand,'U FFT...');
while(first_line_in_block<ifp_params.image_size_pixels(1))
    waitbar((first_line_in_block-1)/ifp_params.image_size_pixels(1),wb_hand);
    line_range = [first_line_in_block ...
        min(ifp_params.image_size_pixels(1),first_line_in_block+max_num_lines-1)];
    if process_in_blocks
        data_block = reader_v_fft.read_chip(line_range, [1 length(new_u)]);
    end
    
    data_block = pfa_fft_zeropad_1d(data_block.', ifp_params.sample_rate);

    writer_final_sicd.write_chip(data_block,[1 first_line_in_block]);
    first_line_in_block=line_range(2)+1;
end
clear writer_final_sicd; % Close file
if process_in_blocks
    reader_v_fft.close(); % Close file
    delete(filename_v_fft); % Intermediate file no longer needed
end

close(wb_hand); % Finally, close waitbar

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////