function ifp_params = select_pulses_samples_cphd( filename, varargin )
%SELECT_PULSES_SAMPLES_CPHD Compute required pulses and samples
% SELECT_PULSES_SAMPLES_CPHD(FILENAME, 'PropertyName', PropertyValue, ...)
% takes a CPHD format phase history file and computers required
% pulse/samples to form an image for a given set of parameters
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
%                         from resolution property.  If resolution is not
%                         givem, then the default is all pulses.
%       sample_range      Set of samples to use for image formation. (For
%                         example, 1:1000.) Default is to calculate this
%                         from resolution property.  If resolution is not
%                         givem, then the default is all samples.
%       channel           Channel to use from CPHD file.  Default is 1.
%       quiet             If false, this reports stats on collection and
%                         IFP parameters.  Default is true.
%
% Output is a structure with the all the input fields completely populated.
%
% Limitations:
% 1) Values for resolution and extent are only approximate.  Does not
% project range vectors into an image plane or consider polar inscription.
% 2) This does not properly calculate required pulses for circular pass
% data.  For any data that spans more than 180 degrees, explicit
% PULSE_RANGE and SAMPLE_RANGE must be given.  Might be nice for circular
% pass datasets if one were allowed to select pulses by center azimuth
% angle (angle from north) and desired resolution.
%
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Handle filenames or reader objects
if ischar(filename)
    ph_reader = open_ph_reader(filename);
elseif isfield(filename,'get_meta') && isfield(filename,'read_cphd')
    ph_reader = filename;
else
    % invalid filename
end
cphd_meta = ph_reader.get_meta();
if ~strcmpi(cphd_meta.CollectionID.RadarMode.ModeType,'SPOTLIGHT')  % The assumptions here are for spotlight data
    error('PFA_FILE:UNSUPPORTED_COLLECT_TYPE','Unsupported collection mode.  Currently only spotlight data is supported.');
end
if strcmpi(cphd_meta.Global.DomainType,'TOA') % Assure spotlight data
    error('SELECT_PULSES_SAMPLES_CPHD:UNSUPPORTED_DOMAIN','TOA domain data not currently supported.');
end

% Extract channel number first
p0 = inputParser;
p0.KeepUnmatched=true;
p0.addParamValue('channel', 1, @(x) isscalar(x)&&x<=cphd_meta.Data.NumCPHDChannels); % Channel from CPHD to process
p0.FunctionName = mfilename;
p0.parse(varargin{:});

% Get the narrowband data.  Note the '[]' in the second parameter of the
% READ_CPHD method means "don't get any sample data" for the pulses.  The
% 'ignore' parameter on the left-hand side is the (nonexistent) pulse data.
[ignore, first_last_pulses_nb] = ph_reader.read_cphd([1 cphd_meta.Data.Channel(p0.Results.channel).NumVectors], [], p0.Results.channel);

% Compute maximum resolutions and scene extents for collect
bw = mean( first_last_pulses_nb.FX2-first_last_pulses_nb.FX1 ); % Total valid bandwidth
bw = mean( first_last_pulses_nb.SCSS*double(cphd_meta.Data.Channel(p0.Results.channel).NumSamples-1));
fc = mean( first_last_pulses_nb.SC0 + ... % Center frequency
    (first_last_pulses_nb.SCSS*double(cphd_meta.Data.Channel(p0.Results.channel).NumSamples)/2) );
freq_step_size = mean(first_last_pulses_nb.SCSS);
LOS        = first_last_pulses_nb.TxPos - first_last_pulses_nb.SRPPos; % Line-of-sight vector between ARP and ORP
[max_resolution, max_extent, delta_azimuth, total_azimuth] = ...
       pulse_info_to_resolution_extent(LOS, fc, freq_step_size, bw, ...
       double(cphd_meta.Data.Channel(p0.Results.channel).NumVectors));
max_resolution = max_resolution * 0.886; % Convert null-to-null to 3dB width

% Parse dependent parameters passed into this function
% We had to wait until here because we wanted default values based on
% what the collection supports. Passed-in values will override the 
% just-computed defaults however.
p1 = inputParser;
p1.KeepUnmatched=true;
p1.addParamValue('resolution',        max_resolution,      @isvector);
p1.addParamValue('pulse_range',       1:cphd_meta.Data.Channel(p0.Results.channel).NumVectors, @isvector);
p1.addParamValue('sample_range',      1:cphd_meta.Data.Channel(p0.Results.channel).NumSamples, @isvector);
p1.addParamValue('quiet',             true,                @islogical); % Display support info.
p1.FunctionName = mfilename;
p1.parse(varargin{:});
ifp_params = p1.Results;
ifp_params.channel = p0.Results.channel;

% If caller specifies a given resolution, rather than a pulse/sample range,
% we'll compute the required pulses/samples and pull both the pulses and
% the samples from the middle of the data.

% Choose samples based on requested range resolution.  We do this first
% since the cross-range resolution depends on the center frequency (we'll
% later use the center frequency of the chosen samples).
c = SPEED_OF_LIGHT; % Speed of light (m/s)
if ~any(strcmp(p1.UsingDefaults,'resolution')) % Resolution specified
    if any(strcmp(p1.UsingDefaults,'sample_range')) % Sample range not specified
        needed_bandwidth        = 0.886*c/(2*ifp_params.resolution(1));
        num_samples             = ceil(needed_bandwidth/freq_step_size);
        ifp_params.sample_range = max(1,floor((cphd_meta.Data.Channel(ifp_params.channel).NumSamples-num_samples)/2)):...
                                  min(cphd_meta.Data.Channel(ifp_params.channel).NumSamples,...
                                  floor((cphd_meta.Data.Channel(ifp_params.channel).NumSamples+num_samples)/2));
    else
        warning('PARSE_IFP_PARAMS:RESOLUTION_IGNORED','Range component of resolution input parameter will be ignored since SAMPLE_RANGE input parameter was passed.')
    end
end
processed_bandwidth = diff(ifp_params.sample_range([1 end]))*freq_step_size;

% Choose pulses base on request cross-range resolution.
processed_fc = mean( first_last_pulses_nb.SC0 + ... % Center frequency of samples to process
    (first_last_pulses_nb.SCSS*sum(ifp_params.sample_range([1 end]))/2) );
if ~any(strcmp(p1.UsingDefaults,'resolution')) % Resolution specified
    if any(strcmp(p1.UsingDefaults,'pulse_range')) % Pulse range not specified
        % Range of pulses to use not given; use all pulses to support the
        % requested resolution (which defaults to "max supported").
        needed_angular_aperture = 0.886*c/(2*processed_fc*ifp_params.resolution(2));
        required_num_pulses_estimate = ceil(needed_angular_aperture/delta_azimuth);
        first_pulse   = max(2,floor((cphd_meta.Data.Channel(ifp_params.channel).NumVectors-required_num_pulses_estimate)/2)+2);
        last_pulse    = min(cphd_meta.Data.Channel(ifp_params.channel).NumVectors-1,first_pulse+required_num_pulses_estimate-2);
        angular_aperture   = 0;
        % Loop to find the minimum number of pulses needed to support the
        % specified resolution.  This gets around having a non-uniform delta
        % theta.
        while (angular_aperture < needed_angular_aperture)&&(first_pulse>1)
            first_pulse = first_pulse - 1;
            last_pulse = last_pulse + 1;
            [ignore, first_last_pulses_nb] = ph_reader.read_cphd([first_pulse last_pulse], [], ifp_params.channel);
            LOS      = first_last_pulses_nb.TxPos - first_last_pulses_nb.SRPPos; % Line-of-sight vector between ARP and ORP
            LOS_norm = repmat(sqrt(sum(LOS.^2,2)),[1 3]);
            LOS_hat  = LOS./LOS_norm;
            angular_aperture = acos(sum(LOS_hat(1,:).*LOS_hat(2,:),2));
        end
        ifp_params.pulse_range = first_pulse:last_pulse;
    else
        warning('PARSE_IFP_PARAMS:RESOLUTION_IGNORED','Cross-range component of resolution input parameter will be ignored since PULSE_RANGE input parameter was passed.')
    end
end

% Given pulse and sample range compute achieved values
[ignore, first_last_pulses_nb] = ph_reader.read_cphd(ifp_params.pulse_range([1 end]), [], ifp_params.channel);
LOS = first_last_pulses_nb.TxPos - first_last_pulses_nb.SRPPos; % Line-of-sight vector between ARP and ORP
[ifp_params.resolution, extent, delta_azimuth, processed_azimuth] = pulse_info_to_resolution_extent(LOS, processed_fc, freq_step_size, processed_bandwidth, length(ifp_params.pulse_range));
ifp_params.resolution = ifp_params.resolution * 0.886; % Convert null-to-null to 3dB width
if ischar(filename)
    ph_reader.close();
end

% Print out some informational data
if ~p1.Results.quiet
    fprintf('%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf('\n');
    if ischar(filename)
        [pathstr, name, ext] = fileparts(filename);
        fprintf('File name "%s%s"\n', name, ext);
        fprintf('\n');
    end
    fprintf('Collection metadata:\n');
    fprintf('   Received bandwidth:           %.1f (MHz)\n', bw/1e6);
    fprintf('   Total collection angle:       %.2f (deg)\n', total_azimuth*180/pi);
    fprintf('   Maximum supported resolutions:\n');
    fprintf('      Range:                     %.3f (m)\n', max_resolution(1));
    fprintf('      Azimuth:                   %.3f (m)\n', max_resolution(2));
    fprintf('   Maximum supported scene extents:\n');
    fprintf('      Range:                     %.0f (m)\n', max_extent(1));
    fprintf('      Azimuth:                   %.0f (m)\n', max_extent(2));
    fprintf('\n');
    fprintf('Processing parameters:\n');
    fprintf('   Range:\n');
    fprintf('      First sample used:         %.0f\n', ifp_params.sample_range(1));
    fprintf('      Last sample used:          %.0f\n', ifp_params.sample_range(end));
    fprintf('      Number of samples:         %.0f\n', length(ifp_params.sample_range));
    fprintf('      Processed RF bandwidth:    %.1f (MHz)\n', processed_bandwidth/1e6);
    fprintf('      Requested resolution:      %.3f (m)\n', p1.Results.resolution(1));
    fprintf('      Achieved resolution:       %.3f (m)\n', ifp_params.resolution(1));
    fprintf('      Extent:                    %.0f (m)\n', max_extent(1));
    fprintf('   Azimuth:\n');
    fprintf('      First pulse used:          %.0f\n', ifp_params.pulse_range(1));
    fprintf('      Last pulse used:           %.0f\n', ifp_params.pulse_range(end));
    fprintf('      Number of pulses:          %.0f\n', length(ifp_params.pulse_range));
    fprintf('      Collection angle:          %.2f (deg)\n', processed_azimuth*180/pi);
    fprintf('      Requested resolution:      %.3f (m)\n', p1.Results.resolution(2));
    fprintf('      Achieved resolution:       %.3f (m)\n', ifp_params.resolution(2));
    fprintf('      Extent:                    %.0f (m)\n', max_extent(2));
    fprintf('\n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%\n');
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////