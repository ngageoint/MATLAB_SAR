function [sicd_meta] = rgazcomp_file(input_filename_cphd, output_filename_sicd, varargin)
%RGAZCOMP_FILE Applies the most simple of IFPs (range-azimuth compression)
%on a file of motion-compensated phase history data in CPHD file format and
%writes to a SICD file with appropriate metadata
%
% RGAZCOMP_FILE(INPUT_FILENAME_CPHD, OUTPUT_FILENAME_SICD, 'PropertyName', PropertyValue, ...)
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
%       sample_rate       Image domain oversample ratios.  Default is 1.
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
% Authors: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters and read metadata
ph_reader = open_ph_reader(input_filename_cphd);
ifp_params = select_pulses_samples_cphd(ph_reader, varargin{:}); % Parse generic IFP parameters
p1 = inputParser;
p1.KeepUnmatched=true;
p1.addParamValue('max_block_size',[], @isscalar);
p1.addParamValue('sample_rate',1, @isscalar);
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
    error('RGAZCOMP_FILE:UNSUPPORTED_COLLECT_TYPE','Unsupported collection mode.  Only spotlight data is supported.');
end
if any(diff(nbdata.SC0)) || any(diff(nbdata.SCSS))
    error('RGAZCOMP_FILE:UNSUPPORTED_COLLECT_TYPE','Unsupported data type.  Must have constant Fx0/Fx_SS.');
end
if any(diff(diff(ifp_params.sample_range)))
    error('RGAZCOMP_FILE:UNSUPPORTED_COLLECT_TYPE','Selected samples must be regularly spaced.');
end
if ~strcmp(cphd_meta.Global.DomainType, 'FX')
    error('RGAZCOMP_FILE:UNSUPPORTED_FILE_TYPE','Unsupported CPHD type.  Currently only FX domain is supported.');
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
    warning('RGAZCOMP_FILE:FILE_PROCESSING','This data must be processed through intermediate files, rather than all in memory.  This could take a while...');
end

ifp_params.image_size_pixels = floor(ifp_params.sample_rate .* [num_samples num_pulses]);
sicd_meta = rgazcomp_sicd_meta(cphd_meta,nbdata,ifp_params);

%% Range FFT across samples
% The FFT will be done in separably in two steps, first in range and then
% in azimuth.

% Setup temporary file
% The range FFT stage includes an in-order read (CPHD is written in
% "pulse-major" order) and an in-order write (data interpolated in range).
% We have to "corner-turn" (or transpose) after the range FFT (as we do on
% all 4 steps), but it is typically faster to read from a file out-of-order
% than write out-of-order, so we will do the corner-turn in the read of the
% azimuth FFT, rather than the write of this range FFT stage.
if process_in_blocks
    filename_range_fft = tempname;
    tempsicdmeta.ImageData.NumRows   = num_pulses;
    tempsicdmeta.ImageData.NumCols   = ifp_params.image_size_pixels(1);
    writer_range_interp = SIOWriter(filename_range_fft, tempsicdmeta, 'include_sicd_metadata', false);
end

% Iterate through sets of pulses for range FFT
first_pulse_in_set=1; % Out of all pulses selected for processing
max_num_lines = floor(max_block_size/(num_samples*data_element_size));
wb_hand=waitbar(0,'Range FFT...');
while(first_pulse_in_set<num_pulses)
    waitbar((first_pulse_in_set-1)/num_pulses,wb_hand);
    pulse_indices = first_pulse_in_set:...
        min(num_pulses,first_pulse_in_set+max_num_lines-1);
    data_block = ph_reader.read_cphd(ifp_params.pulse_range(pulse_indices),...
        ifp_params.sample_range, ifp_params.channel);
    
    data_block = pfa_fft_zeropad_1d(data_block,ifp_params.sample_rate);

    if process_in_blocks
        writer_range_interp.write_chip(data_block,[1 first_pulse_in_set]); % Save range interpolated data in this block
    end
    first_pulse_in_set=pulse_indices(end)+1;
end
if process_in_blocks
    clear writer_range_interp % Close file
end
ph_reader.close();

%% Azimuth FFT across vectors
if process_in_blocks
    reader_range_fft = open_reader(filename_range_fft);
end
writer_final_sicd = SICDWriter(output_filename_sicd, sicd_meta); % Final image will be written here

% Iterate through blocks of data for azimuth FFT
first_line_in_block=1; % Out of all rows in the final image
max_num_lines = floor(max_block_size/(num_pulses*data_element_size));
waitbar(0,wb_hand,'Azimuth FFT...');
while(first_line_in_block<ifp_params.image_size_pixels(1))
    waitbar((first_line_in_block-1)/ifp_params.image_size_pixels(1),wb_hand);
    line_range = [first_line_in_block ...
        min(ifp_params.image_size_pixels(1),first_line_in_block+max_num_lines-1)];
    if process_in_blocks
        data_block = reader_range_fft.read_chip(line_range, [1 num_pulses]);
    end
    
    data_block = pfa_fft_zeropad_1d(data_block.',ifp_params.sample_rate);
    if isfield(sicd_meta,'SCPCOA') && isfield(sicd_meta.SCPCOA,'SideOfTrack') && ...
            isequal(sicd_meta.SCPCOA.SideOfTrack, 'R')
        data_block = flipud(data_block); % View from above
    end

    writer_final_sicd.write_chip(data_block,[1 first_line_in_block]);
    first_line_in_block=line_range(2)+1;
end
clear writer_final_sicd; % Close file
if process_in_blocks
    reader_range_fft.close(); % Close file
    delete(filename_range_fft); % Intermediate file no longer needed
end

close(wb_hand); % Finally, close waitbar

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////