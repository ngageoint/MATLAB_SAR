function bp_file( infilename, outfilename, varargin )
%BP_FILE Forms an image from CPHD data via back projection 
% BP_FILE(infilename, outfilename, 'PropertyName', PropertyValue, ...)
% takes a CPHD format phase history file (infilename) and forms an image,
% writing to the SICD complex image file specified (outfilename).
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
%       sample_rate       Samples per IPR.  Default is 1.5.
%       channel           Channel to use from CPHD file.  Default is 1.
%       image_size_meters Scene size for image formation ([range_size
%                         azimuth_size]) in meters.  Default is total size
%                         supported by collect.
%       center            Center of image to form in ECF coordinates.
%                         Default is scene reference point.
%       grid_type         'slant', 'ground', or 'DEM'.  Default is ground.
%       max_image_block_size Largest size of subimage to be processed in
%                         memory at any given time (in bytes)
%       max_pulse_block_size Maximum amount of phase history data to be
%                         held in memory at any given time (in bytes).
%       quiet             If false, this reports stats on collection and
%                         IFP parameters.  Default is true.
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
% Input parameters specific to bp_file.m
if ispc
    [uv, sv] = memory;
    % max_*_block_size is just the size of the largest single array we
    % will handle at once.  We will need significantly more memory than
    % this.  A factor of 64 seems to be reasonable.  Its a nice round
    % power of two that seems to avoid swapping in limited testing.
    default_max_block_size = min(sv.PhysicalMemory.Available/64, ...
        uv.MaxPossibleArrayBytes);
else % Can't gauge memory in UNIX, so just put arbitrary value
    default_max_block_size = 2^27; % About 100Meg
end
p1 = inputParser;
p1.KeepUnmatched=true;
p1.addParamValue('max_image_block_size', default_max_block_size, @isscalar);
p1.addParamValue('max_pulse_block_size', default_max_block_size, @isscalar);
p1.addParamValue('show_waitbar', true, @isscalar);
p1.FunctionName = mfilename;
p1.parse(varargin{:});
ph_reader = open_ph_reader(infilename);
% Select pulses with which to form image
ifp_params = setstructfields(p1.Results, select_pulses_samples_cphd(ph_reader, varargin{:}));
% Read in per-pulse metadata (with no sample data yet)
[ignore, nbdata] = ph_reader.read_cphd(ifp_params.pulse_range, [], ifp_params.channel);
% Compute grid of points onto which we will form the image
num_pulses_used  = length(ifp_params.pulse_range);
num_samples_used = length(ifp_params.sample_range);
ifp_params = setstructfields(ifp_params, bp_parse_grid_params(nbdata, num_samples_used, varargin{:}));

% We're assuming neither the entire output image nor the entire set of
% pulses can fit in memory.  Here we get the size of sub-images and the
% number of pulses in a pulse "block". The 'min' below lets the entire
% image or set of pulses reside in memory if they're smaller than the
% defined max{Image|Pulse}BlockSize (defined above).
% For now, we process the image in sets of full rows, since that's the
% order that image is written in.  This allows us to write in order to
% file.
imageBlockSize   = ifp_params.image_size_pixels;
bytes_per_row    = ifp_params.image_size_pixels(2)*4; % Because we have single precision image data
if isempty(ifp_params.grid) % Can't currently iterate through chips of arbitrary grid
    imageBlockSize(1)= min( imageBlockSize(1), floor(ifp_params.max_image_block_size/bytes_per_row) );
end
% Number of pulses that will be read into memory at any given time
pulseBlockSize   = min( num_pulses_used, floor(ifp_params.max_pulse_block_size/...
    (num_samples_used*4)) ); % Assume 4 byte samples, largest allowed by CPHD
% This is the number of sub-images and pulse "blocks" we'll be using.  Note
% the 'numImageBlocks' is a 2-vector indicating the number of sub images in
% the row and column direction.  The total number of sub-images we'll be
% forming is the numImageBlocks(1)*numImageBlocks(2).  'numPulseBlocks' is
% just a scalar value.
numImageBlocks = ceil(ifp_params.image_size_pixels./imageBlockSize);
numPulseBlocks = ceil(num_pulses_used./pulseBlockSize);
totalNumBlocks = prod(numImageBlocks)*numPulseBlocks;

%% Form the image in chips using sets of pulses.
% First setup image file to which we will write
sicdmeta = bp_sicd_meta(ph_reader.get_meta(),nbdata,ifp_params); % Form SICD metadata
imageWriter = SIOWriter(outfilename, sicdmeta); % Open output file

% Setup the GUI waitbar (if enabled).  Note, we'll set it up in one of two
% places.  If we're doing a single image chip and single pulse block we'll
% let 'bpBasic' handle it, otherwise we'll do it here at the higher level
% (since we're doing multiple calls to 'bpBasic').
showWaitbarHere = ifp_params.show_waitbar  &&  (max(numImageBlocks) > 1  ||  numPulseBlocks > 1);
if showWaitbarHere
    wb = waitbar(0); tic
end

% Loop over image chips forming the complete image chip at the desired
% resolution and writing it to the correct location in the output image
% file.  We'll build the chip by forming "sub-images" from sets of pulses
% (that fit in memory) and accumulating the resulting image.  When we've
% processed all the desired pulses we'll have a complete image chip.
blockCount = 0;
for imageBlockRow = 1:numImageBlocks(1)
    for imageBlockCol = 1:numImageBlocks(2) % This loop should always be a single iteration, since we process sets of full rows.
        imageBlockStart = [imageBlockRow-1,imageBlockCol-1] .* imageBlockSize + [1,1];
        imageBlockEnd   =  min( ifp_params.image_size_pixels, imageBlockStart + imageBlockSize - [1,1]);
        if ~isempty(ifp_params.grid) % Arbitrary grid.  Should be only single image block.  Just iterating over pulses.
            data.x_mat = ifp_params.grid(:,:,1);
            data.y_mat = ifp_params.grid(:,:,2);
            data.z_mat = ifp_params.grid(:,:,3);
        elseif ~strcmpi(ifp_params.grid_type, 'DEM') % Plane
            [x_grid,y_grid] = ndgrid( ifp_params.col_coords(imageBlockStart(2):imageBlockEnd(2)), ...
                ifp_params.row_coords(imageBlockStart(1):imageBlockEnd(1)) );
            data.x_mat = ifp_params.center(1) + ifp_params.col_unit_vector(1)*x_grid + ifp_params.row_unit_vector(1)*y_grid;
            data.y_mat = ifp_params.center(2) + ifp_params.col_unit_vector(2)*x_grid + ifp_params.row_unit_vector(2)*y_grid;
            data.z_mat = ifp_params.center(3) + ifp_params.col_unit_vector(3)*x_grid + ifp_params.row_unit_vector(3)*y_grid;
        else % DEM,  Form lat/long grid with heights extracted from DEM.
            lat_subset = ifp_params.row_coords(imageBlockStart(1):imageBlockEnd(1));
            lon_subset = ifp_params.col_coords(imageBlockStart(2):imageBlockEnd(2));
            block_size = [length(lon_subset) length(lat_subset)];
            [data.x_mat, data.y_mat, data.z_mat]=deal(zeros(block_size));
            for i=1:block_size(1)
                for j=1:block_size(2)
                    h = ifp_params.F_1(lat_subset(j), lon_subset(i));
                    pos_lla = [lat_subset(j), lon_subset(i), h];
                    pos_ecef = geodetic_to_ecf(pos_lla);
                    data.x_mat(i,j) = pos_ecef(1);
                    data.y_mat(i,j) = pos_ecef(2);
                    data.z_mat(i,j) = pos_ecef(3);
                end
            end
        end
        
        image_chip = zeros(size(data.x_mat));
        for pulseBlock = 1:numPulseBlocks
            blockCount = blockCount + 1;
            if showWaitbarHere
                t_sofar = toc;
                t_est = (t_sofar*(totalNumBlocks+1)/blockCount)-t_sofar;
                if blockCount==1
                    tr_message = 'calculating...';
                else
                    tr_message = datestr(datenum(0,0,0,0,0,t_est),13);
                end
                wb_message=sprintf('Processing block %d of %d, Time remaining: %s', ...
                    blockCount, totalNumBlocks, tr_message);
                waitbar((blockCount-1)/(totalNumBlocks),wb,wb_message);
            end
            
            first_pulse = (pulseBlock-1)*pulseBlockSize + 1;
            last_pulse  = min(first_pulse + pulseBlockSize - 1, num_pulses_used);
            [phase_history, nbdata] = ph_reader.read_cphd(...
                ifp_params.pulse_range(first_pulse:last_pulse), ifp_params.sample_range, ifp_params.channel);
            
            %% Call the backprojection function with the appropriate inputs
            % Convert metadata into structure that bpBasic function wants
            % Antenna positions x, y, and z (all in meters)
            data.Tx.X    = nbdata.TxPos(:,1,1).';
            data.Tx.Y    = nbdata.TxPos(:,2,1).';
            data.Tx.Z    = nbdata.TxPos(:,3,1).';
            data.Rcv.X    = nbdata.RcvPos(:,1,1).';
            data.Rcv.Y    = nbdata.RcvPos(:,2,1).';
            data.Rcv.Z    = nbdata.RcvPos(:,3,1).';
            data.R0      = (sqrt(sum((nbdata.TxPos-nbdata.SRPPos).^2,2)) + ... % Range from antenna to motion comp point
                           sqrt(sum((nbdata.RcvPos-nbdata.SRPPos).^2,2)))/2;
            data.minF    = nbdata.SC0 + (nbdata.SCSS * double(ifp_params.sample_range(1))); % Minimum frequency for each pulse (Hz)
            data.deltaF  = nbdata.SCSS; % Frequency step size (Hz)
            data.phdata  = phase_history;
            data.hide_waitbar = showWaitbarHere;
            % The following is the only line of actual IFP in this whole function:
            complex_image = bpBasic(data);
            % "Basebanding".  This centers the frequency support, even if
            % the image is formed to a DEM.  When working on image domain
            % data, we need an average center frequency across all pulses.
            f_c          = mean(nbdata.SC0  + (nbdata.SCSS * ...
                double(ifp_params.sample_range(ceil(length(ifp_params.sample_range)/2))-1)));
            c = SPEED_OF_LIGHT; % Speed of light (m/s)
            complex_image = complex_image.*...
                exp(-1i*4*pi*f_c/c*data.x_mat*ifp_params.range_vect_hat(1)).*...
                exp(-1i*4*pi*f_c/c*data.y_mat*ifp_params.range_vect_hat(2)).*...
                exp(-1i*4*pi*f_c/c*data.z_mat*ifp_params.range_vect_hat(3));
            
            image_chip = image_chip + complex_image;
        end
        imageWriter.write_chip(image_chip,imageBlockStart([2 1]));
    end
end

%% Clean up
if showWaitbarHere
    close(wb);
end
ph_reader.close();
clear imageWriter % Close output image file

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////