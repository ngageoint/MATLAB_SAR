function [image_information_metric, rniirs, diagnostic] = RGIQE(filename, varargin)
% RGIQE Radar generalized image quality equation
%    [image_information_metric, rniirs] = RGIQE(bandwidth, nesz_db, graze_degrees)
%    [image_information_metric, rniirs] = RGIQE(filename, 'PropertyName', PropertyValue, ...)
%
% An implementation of a generalized image quality equation to provide a
% quantitative image quality metric and an RNIIRS prediction.
%
% RNIIRS (Radar National Imagery Interpretation Rating Scale) is an image
% quality scale from 1 to 9.  This "quality" refers to local image fidelity
% and a users ability to perform certain detection and recognition tasks
% with an image.  The "quality" that RNIIRS attempts to capture does NOT
% refer to image size, geopositioning accuracy, radiometric accuracy, or
% any other host of characterists that may improve an image's utility.
% Although RNIIRS is subjective, generally assessed by a trained analyst,
% an image quality equation (IQE) may be used to predict RNIIRS for a given
% SAR collect.
%
% The General Image Quality Equation (GIQE) provided here uses generic,
% fundamental, sensor-indepenedent parameters of radar collection to
% measure local image fidelity by quantifying the information density in
% the ground in terms of bits\m^2, using Shannon's channel capacity
% theorem.  Unlike RNIIRS, this value is quantitative and not subjective or
% based on analysts' assessments.  We can however compute an RNIIRS
% prediction based off of it by fitting the information metric to analysts
% ratings.  Both the information theory metric and the RNIIRS prediction
% are computed in this function.
%
%    [image_information_metric, rniirs] = RGIQE(bandwidth, nesz_db, graze_degrees)
% INPUTS:
%
%    bandwidth               Transmitted (and processed) bandwidth in Hz.
%    nesz_db                 Noise equivalent sigma-0 in dB.  Should
%                               capture all noise, including multiplicative
%                               sources.
%    graze_degrees           Grazing angle in degrees.
%
%
%    [image_information_metric, rniirs] = RGIQE(filename, 'PropertyName', PropertyValue, ...)
% INPUTS:
%
%    filename                Must either be a complex image in a format
%                               recognized by the MATLAB SAR Toolbox or an
%                               XML file with SICD-formatted XML.  If empty
%                               or not passed, the user will be prompted to
%                               select a file with a dialog box.
%
%    Property name           Description
%    multi_noise             Multiplicative noise.  This value is equal to:
%
%                               MNR * sigma_0_background
%
%                               where MNR is the multiplicative noise ratio
%                               and sigma_0_background is the assumed power
%                               level of the background and ambiguous
%                               areas.
%                               MNR can further be broken down into:
%
%                               MNR = ISLR + QNR + AMBR
%
%                               where ISLR is the integrated sidelobe
%                               ratio, QNR is the quantization noise ratio,
%                               and AMBR is the ambiguity ratio.
%                               (Default = 0)
%    interactive             "AUTO", "NOISE", "SIGNAL", "BOTH", or "NONE".
%                               Enable interactive measurement of noise and
%                               signal levels.  "AUTO" only runs
%                               interactive measurment if metadata from
%                               source file does not provide it. (Default
%                               is "AUTO").
%    frames                  If a file contains more than one complex
%                               dataset, this value specifies for which of
%                               those datasets to compute metrics. (Default
%                               is all datasets in the file.)
%    signal_sigma_override   Allows caller to override signal level used
%                               for computation (and INTERACTIVE input
%                               parameter).  In general, this is not
%                               recommended, as RNIIRS estimation function
%                               was trained using the default values. If
%                               the dataset provided is radiometrically
%                               calibrated, this value will be treated as
%                               calibrated signal level in the ground plane
%                               (sigma-0).  If the dataset is not
%                               calibrated, this value will be assume to be
%                               average pixel power. (Default is 1 square
%                               meter of radar return per square meter on
%                               the ground for co-pol, 0.25 for cross-pol.)
%
% Note that the selection of signal level and noise level are fundamentally
% different in this model.  Noise is measureable, quantifiable, and likely
% fairy accuractely predictable (at least for additive noise).  There is a
% right answer for what the noise level is for any given dataset. On the
% other hand, signal level is somewhat arbitrarily chosen, selected to be
% the power level of the scene content that one may be most interested in
% examining.  It is for this reason that this function allows for a signal
% level override, but not a direct noise level override (although noise can
% be adjusted with the MULTI_NOISE input parameter.)
%
% Author: Wade Schwartzkopf, NGA/Research

% Of the lines in this file, >95% is file and input parameter parsing and
% documention, and <5% actually computing the image quality metrics...

%% Bandwidth formulation
% Assumes azimuth bandwidth to match range in ground plane
% NESZ provided should capture all noise, including multiplicative sources
if nargin == 3 && isnumeric(filename) && isnumeric(varargin{1}) && isnumeric(varargin{2})
    bandwidth = filename;  % Hz
    nesz = varargin{1};  % in dB
    graze = varargin{2};  % in degrees
    bandwidth = bandwidth*2/SPEED_OF_LIGHT;  % Convert units from Hz to cycles/m
    bandwidth = bandwidth*cosd(graze);  % Convert from slant plane to ground plane
    bandwidth = bandwidth.^2;  % Convert from 1D bandwidth to 2D bandwidth area (assumes azimuth bandwidth same as range)
    image_information_metric = bandwidth.*log2(1+1/(10^(nesz/10)));
    rniirs = estimate_rniirs(image_information_metric);
    diagnostic = struct();  % No diagnostics for this version
    return;
end

%% SICD/SLC version.
%% Open input file
if ((nargin<1)||isempty(filename)) % If no filename was give, use dialog box to ask for one
    % Recall last interactively selected path used
    if ispref('matlab_sar_toolbox','last_used_directory')
        pathname = getpref('matlab_sar_toolbox','last_used_directory');
        if ~ischar(pathname)||~exist(pathname,'dir')
            pathname = pwd;
        end
    else
        pathname = pwd;
    end
    [filename,pathname]=uigetfile(sar_file_extensions(),...
        'Open SAR Data File',pathname);
    setpref('matlab_sar_toolbox','last_used_directory',pathname);
else % Path was already passed in with filenames
    pathname='';
end
if(filename)
    fullfilename=[pathname filename];
else % filename=0.  Cancel was pressed, instead of a file being chosen.
    return;
end
if ~isempty(pathname)
    setpref('matlab_sar_toolbox','last_used_directory',pathname); %store path
end
try
    ro = open_reader(fullfilename);
    % Check for multiple images
    if ~iscell(ro)
        ro = {ro};
    end
    meta = cell(numel(ro),1);
    for i = 1:numel(ro)
        meta{i} = ro{i}.get_meta();
    end
catch
    try
        fid = fopen(fullfilename);
        meta = {sicdxml2struct(xmlread( java.io.StringBufferInputStream( char(fread(fid)'))))};
        fclose(fid);
    catch
        error('RGIQE:InvalidFile','Invalid File.');
    end
end

%% Parse other input parameters
p = inputParser;
p.addParameter('multi_noise', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('interactive', 'AUTO', @(x) ismember(upper(x), ...
    {'AUTO', 'NOISE', 'SIGNAL', 'BOTH', 'NONE'}));
p.addParameter('frames', 1:numel(meta), @(x) all(ismember(x, 1:numel(meta))));
p.addParameter('signal_sigma_override', [], @(x) isnumeric(x) && isscalar(x));
parse(p,varargin{:})

%% Iterate through datasets provided
meta = meta(p.Results.frames);
image_information_metric = zeros(numel(meta),1);
rniirs = zeros(numel(meta),1);
for i = 1:numel(meta)
    meta{i} = derived_sicd_fields(meta{i});  % Assure all radiometric fields are derived if not there

    %% Get noise level
    if isfield(meta{i}, 'Radiometric') && ...
            isfield(meta{i}.Radiometric,'NoiseLevel') && ...
            isfield(meta{i}.Radiometric.NoiseLevel,'NoiseLevelType') && ...
            strcmpi(meta{i}.Radiometric.NoiseLevel.NoiseLevelType,'ABSOLUTE') && ...
            isfield(meta{i}.Radiometric.NoiseLevel,'NoisePoly') && ...
            ~ismember(upper(p.Results.interactive), {'NOISE', 'BOTH'})
        % SICD metadata can only tell us additive noise.  Estimation of
        % multiplicative noise must come from elsewhere.
        noise_sigma = 10^(meta{i}.Radiometric.NoiseLevel.NoisePoly(1)/10) + p.Results.multi_noise;
    elseif ismember(upper(p.Results.interactive), {'AUTO', 'NOISE', 'BOTH'})
        % Measured noise contains both additive and multiplicative noise
        uiwait(msgbox('Please select a noise region for measurement.','Noise measurement'));
        aoi = mitm_viewer(fullfilename,'mode','aoi','closeAfterSelect',true,'initialFrame',i);
        if isempty(aoi), return; end  % Cancel was selected
        cdata = double(ro{i}.read_chip(...
            [aoi(1) aoi(1)+aoi(3)-1],[aoi(2) aoi(2)+aoi(4)-1]));
        noise_sigma = mean(abs(cdata(:)).^2);
    else
        noise_sigma = NaN;
    end
    if isfield(meta{i}, 'Radiometric') && ...
            isfield(meta{i}.Radiometric, 'SigmaZeroSFPoly')
        % Scale to be equal to noise equivalent sigma zero
        noise_sigma = meta{i}.Radiometric.SigmaZeroSFPoly(1) * noise_sigma;
    end

    %% Get signal level
    % One issue with using the Shannon-Hartley channel capacity formulation
    % is that it requires a signal level.  Certainly the amount of
    % information a SAR image conveys is related to its signal level, so
    % this is how we hope it would behave.  A target with brighter returns
    % is always easier to interpret than a target that is in the noise.
    % However, it presents a problem to image quality prediction, since the
    % signal level can't be known prior to measurement.  This constraint is
    % implicitly true for traditional RNIIRS assessments as well, which
    % assumes specific types of targets.
    % One approach for setting the signal level used in the channel
    % capacity equation would be to assume a fixed radiometric signal level
    % that is representative of the sorts of scene content we are generally
    % interested in and use that same signal level across all images for
    % which we compute this metric.  This allows for an even
    % apples-to-apples comparison across images of varying content and is
    % consistent with the traditional RNIIRS viewpoint, which assumes
    % specific types of tasks and scene content.
    % Another approach would be to use the signal level in the actual scene
    % by measuring the average pixel power across a full image or AOI. This
    % approach may be more accurate at measuring the true information in
    % any given image area.  However, it may be less consistent with the
    % traditional RNIIRS viewpoint, and this approach can only be used
    % after the collection has been made, not prior to it for planning
    % collection and predicting image quality.
    % We also allow the caller of the function to pass their own signal
    % level. Passing a non-default signal level is not recommended for
    % RNIIRS estimation, since the RNIIRS function was trained on the
    % default values.  However, we provide this flexibility in case it
    % might be educational to see how it affects the information theory
    % metric.
    if ~ismember('signal_sigma_override', p.UsingDefaults)
        signal_sigma = p.Results.signal_sigma_override;
    elseif isfield(meta{i}, 'Radiometric') && ...
            isfield(meta{i}.Radiometric, 'SigmaZeroSFPoly') && ...
            ~ismember(upper(p.Results.interactive), {'SIGNAL', 'BOTH'})
        % Default value: 1 square meter of return per square meter on
        % the ground (roughly similar to a large vehicle over its area).
        % Ulaby plots show cross-pol generally 5-10dB lower than
        % co-pol.  We pick 1/4 as simple round number.
        if isfield(meta{i},'ImageFormation') && ...
                isfield(meta{i}.ImageFormation,'TxRcvPolarizationProc') && ...
                numel(unique(split(meta{i}.ImageFormation.TxRcvPolarizationProc,':')))>1
            signal_sigma = 0.25;
        else
            signal_sigma = 1;
        end
    elseif ismember(upper(p.Results.interactive), {'AUTO', 'SIGNAL', 'BOTH'})
        ButtonName = questdlg(['Do you want to use average signal '...
            'power across entire image or select an AOI?'], ...
            'Signal level estimation', ...
            'Image', 'AOI', 'Image');
        if isempty(ButtonName), return; end %window was closed CANCEL
        if strcmpi(ButtonName,'AOI')
            aoi = mitm_viewer(fullfilename,'mode','aoi','closeAfterSelect',true,'initialFrame',i);
            if isempty(aoi), return; end  % Cancel was selected
            cdata = double(ro{i}.read_chip(...
                [aoi(1) aoi(1)+aoi(3)-1],[aoi(2) aoi(2)+aoi(4)-1]));
            signal_sigma = mean(abs(cdata(:)).^2) - noise_sigma;
        else
            samplesize=[1000 1000]; % Exract 1000x1000 array of samples to estimate mean
            subsample=ceil(double([meta{i}.ImageData.NumCols meta{i}.ImageData.NumRows])./samplesize);
            data = abs(single(ro{i}.read_chip([1 meta{i}.ImageData.NumCols],...
                    [1 meta{i}.ImageData.NumRows],subsample)));
            signal_sigma = mean(abs(data(:)).^2) - noise_sigma;
        end
        if isfield(meta{i}, 'Radiometric') && ...
                isfield(meta{i}.Radiometric, 'SigmaZeroSFPoly')
            % Scale to be equal to noise equivalent sigma zero
            signal_sigma = meta{i}.Radiometric.SigmaZeroSFPoly * signal_sigma;
        end
        warning('RGIQE:SignalLevelMeasured',['More uncertainty with ' ...
            'RNIIRS estimates exists when signal level is sampled, ' ...
            'rather than when using a radiometric constant.']);
    else
        signal_sigma = NaN;
    end

    %% Save some diagnostic statistics
    if isfield(meta{i}, 'Radiometric') && ...
            isfield(meta{i}.Radiometric, 'SigmaZeroSFPoly')
        diagnostic.nesz(i) = 10*log10(noise_sigma);
    end
    diagnostic.snr(i) = signal_sigma/noise_sigma;
    diagnostic.sp_resolution(i) = sqrt(meta{i}.Grid.Row.ImpRespWid * ...
        meta{i}.Grid.Col.ImpRespWid);
    diagnostic.ellipicity(i) = max(meta{i}.Grid.Row.ImpRespWid, ...
        meta{i}.Grid.Col.ImpRespWid)/min(meta{i}.Grid.Row.ImpRespWid, ...
        meta{i}.Grid.Col.ImpRespWid);
    if numel(meta)>1
        diagnostic.meta{i} = meta{i};
    else
        diagnostic.meta = meta{i};
    end

    %% Compute bandwidth area
    % The cosine of the slope angle is proper scale factor to project the
    % bandwidth area into the ground plane.
    diagnostic.bandwidth_area(i) = ...
        meta{i}.Grid.Col.ImpRespBW * meta{i}.Grid.Row.ImpRespBW * ...
        cosd(meta{i}.SCPCOA.SlopeAng);
    % Computing area from vertices is much more work, but it may help
    % visualize.
    % Project slant plane spatial frequency bounds to ground plane
    % coords_slant_x = meta{i}.Grid.Col.KCtr + ...
    %     meta{i}.Grid.Col.ImpRespBW * [-1 1 1 -1 -1]/2;
    % coords_slant_y = meta{i}.Grid.Row.KCtr + ...
    %     meta{i}.Grid.Row.ImpRespBW * [1 1 -1 -1 1]/2;
    % ruvect = [meta{i}.Grid.Row.UVectECF.X meta{i}.Grid.Row.UVectECF.Y meta{i}.Grid.Row.UVectECF.Z];
    % cuvect = [meta{i}.Grid.Col.UVectECF.X meta{i}.Grid.Col.UVectECF.Y meta{i}.Grid.Col.UVectECF.Z];
    % coords_slant_3d = ruvect'*coords_slant_y + cuvect'*coords_slant_x;
    % Compute basis vectors for ground plane
    % gpn = wgs_84_norm([meta{i}.GeoData.SCP.ECF.X meta{i}.GeoData.SCP.ECF.Y meta{i}.GeoData.SCP.ECF.Z])';
    % gruvect = project_ground(ruvect, gpn); % Project range vector to ground plane
    % gruvect = gruvect/norm(gruvect);
    % gcuvect = cross(gruvect,gpn);
    % gcuvect = gcuvect/norm(gcuvect);
    % Project spatial frequencies to ground plane
    % coords_ground_3d = coords_slant_3d;
    % coords_ground_x = zeros(size(coords_slant_x));
    % coords_ground_y = zeros(size(coords_slant_y));
    % for j=1:numel(coords_slant_x)
    %     % Planes in ECF
    %     coords_ground_3d(:,j) = project_ground(coords_slant_3d(:,j)',gpn)';
    %     % Convert to 2D basis vectors in ground plane
    %     coords_ground_x(j) = dot(coords_ground_3d(:,j),gcuvect);
    %     coords_ground_y(j) = dot(coords_ground_3d(:,j),gruvect);
    % end
    % This is a way to compute the area of a any polygon from its vertices.
    % diagnostic.bandwidth_area = sum((coords_ground_x(1:(end-1)).*coords_ground_y(2:end)) - ...
    %     (coords_ground_y(1:(end-1)).*coords_ground_x(2:end)))/2;
    % Visualize spatial frequency projections in slant and ground planes
    % figure;
    % title('Spatial frequency bounds 3D');
    % g_plane = gruvect'*[100 100 -100 -100 100] + gcuvect'*[-100 100 100 -100 -100];
    % plot3(coords_slant_3d(1,:),coords_slant_3d(2,:),coords_slant_3d(3,:),...
    %     coords_ground_3d(1,:),coords_ground_3d(2,:),coords_ground_3d(3,:),...
    %     g_plane(1,:),g_plane(2,:),g_plane(3,:));
    % legend('Slant plane bounds', 'Ground plane bounds', 'Ground plane')
    % figure;
    % plot(coords_slant_x, coords_slant_y, coords_ground_x, coords_ground_y);
    % title('Spatial frequency bounds');
    % legend('Slant plane bounds', 'Ground plane bounds')

    %% After 300+ line of setup code and comments, we finally compute quality metric...
    % Shannon-Hartley theorem of channel capacity
    % C = B * log2(1 + (S/N))
    % We propose a 2-dimension information rate based on a 2-dimensional
    % bandwidth.  This 2D bandwidth is an area in a space where each axis
    % is cycles/m, as opposed to 1D distance in Hz (cycles per second).
    % Therefore the units of this metric are bits/m^2.
    image_information_metric(i) = diagnostic.bandwidth_area(i) * log2(1 + diagnostic.snr(i));

    %% Map bits/m^2 into RNIIRS
    if ismember('signal_sigma_override', p.UsingDefaults)
        rniirs(i) = estimate_rniirs(image_information_metric(i));
    else  % If assumed signal is not default, this RNIIRS mapping is not valid
        rniirs(i) = NaN;
    end
end
try
    for i = 1:numel(meta)
        ro{i}.close();
    end
end

end

function rniirs = estimate_rniirs(image_information_metric)
    % The coefficients here were derived by fitting hundreds of SAR
    % datasets with varying parameters to analysts ratings.  This
    % author prefers the information metric but also provides this as
    % well for comparison to this legacy scale.
    coeffs = [.3960 3.7555];  % Probably more digits of precision than we really have
    rniirs = polyval(coeffs,log2(image_information_metric));
    % For completeness, we handle the very low RNIIRS case where the
    % logarithmic equation would have been negative.  Negative is not
    % allowed in the RNIIRS scale. We compute where a line tangent to
    % the bits/m^2-to-RNIIRS function intersects with the origin.  We
    % will use this tangent line as the mapping for very low RNIIRS,
    % since it approaches zero as information approaches zero and
    % remains non-negative. This way the entire mapping function
    % between bits/m^2 and RNIIRS is continuous and has a continuous
    % derivative.  The function will be logarithmic with respect to the
    % information metric above this point (presumably nearly all data
    % will fit in here) and linear with respect to the metric below
    % this point (extremely low RNIIRS).  Note that since NGA has no
    % precedent or analysts ratings for RNIIRS<1, we are really free to
    % define this however we want.
    linlog_transition = exp(1 - coeffs(2)*log(2)/coeffs(1));
    % Derivative of linear portion (as well as intersection with log portion)
    d_rniirs = coeffs(1)/(log(2)*linlog_transition);
    % This is what the two lines would look line at very low RNIIRS:
    % x=linspace(0,.02,1000);
    % figure; plot(x, polyval(coeffs,log2(x)), x, x*d_rniirs);
    % That was a lot of explanation for a case that will likely rarely
    % be used, wasn't it?
    if image_information_metric<linlog_transition
        rniirs = d_rniirs * image_information_metric;
    end
end

function out = project_ground(point, gpn)
    % Projection of a point along a given direction to a plane is
    % just the intersection of the line defined by that point (l0)
    % and direction (l) and the plane defined by a point in the
    % plane (p0) and the normal (p):
    % l0 - ((l0 - p0).p/(l.p))*l
    % where . represents the dot product.
    % In this specific case, the origin is in the plane, and the projection
    % direction is the same as the plane normal, so this equation is much
    % simplified:
    out = point - dot(point,gpn)*gpn;
end