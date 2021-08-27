function [complex_image, sicdmeta, grid] = bp_mem(filename, varargin )
%BP_MEM Forms an image from supplied phase history data via the back projection algorithm 
% BP_MEM(FILENAME, 'PropertyName', PropertyValue, ...) takes a CPHD
% format phase history file and forms an image with all processing done in
% memory.
%
% Assumes all phase history and output complex data can fit in memory.
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
%       sample_rate       Samples per IPR.  Default is 1.5.
%       pulse_range       Set of pulses to use for image formation. (For
%                         example, 1:1000.) Default is to calculate this
%                         from resolution property.
%       sample_range      Set of samples to use for image formation. (For
%                         example, 1:1000.) Default is to calculate this
%                         from resolution property.
%       channel           Channel to use from CPHD file.  Default is 1.
%       quiet             If false, this reports stats on collection and
%                         IFP parameters.  Default is true.
%       center            Center of image to form in ECF coordinates.
%                         Default is scene reference point.
%       image_size_meters Scene size for image formation ([range_size
%                         azimuth_size]) in meters.  Default is total size
%                         supported by collect.
%       grid              The grid of image points (ECEF) onto which the
%                         image should be formed.  The format of the grid
%                         parameters is identical to the (MxNx3) grid
%                         returned from this routine if no grid is
%                         supplied.
%       grid_type         'slant', 'ground', or 'DEM'.  Default is ground.
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Initial setup of useful values
% Parse input parameters
ph_reader = open_ph_reader(filename);
ifp_params = select_pulses_samples_cphd(ph_reader, varargin{:});
% Read in the data.  The phase history data will have pulses in the columns
% (i.e. each column of the returned PH is a pulse).
cphd_meta = ph_reader.get_meta();
[phase_history, nbdata] = ph_reader.read_cphd(ifp_params.pulse_range, ifp_params.sample_range, ifp_params.channel);
ph_reader.close();

%% Setup imaging grid if none was supplied
grid_params = bp_parse_grid_params( nbdata, length(ifp_params.sample_range), varargin{:} );
if isempty(grid_params.grid)
    if ~strcmpi(grid_params.grid_type, 'DEM') % Plane
        [x_grid,y_grid] = ndgrid(grid_params.col_coords,grid_params.row_coords);
        data.x_mat=grid_params.center(1) + ...
            grid_params.col_unit_vector(1)*x_grid + grid_params.row_unit_vector(1)*y_grid;
        data.y_mat=grid_params.center(2) + ...
            grid_params.col_unit_vector(2)*x_grid + grid_params.row_unit_vector(2)*y_grid;
        data.z_mat=grid_params.center(3) + ...
            grid_params.col_unit_vector(3)*x_grid + grid_params.row_unit_vector(3)*y_grid;
    else % DEM
        [data.x_mat, data.y_mat, data.z_mat]=deal(zeros(grid_params.image_size_pixels([2 1])));
        for i=1:grid_params.image_size_pixels(2)
            for j=1:grid_params.image_size_pixels(1)
                h = grid_params.F_1(grid_params.row_coords(j), grid_params.col_coords(i));
                pos_lla = [grid_params.row_coords(j), grid_params.col_coords(i), h];
                pos_ecef = geodetic_to_ecf(pos_lla);
                data.x_mat(i,j) = pos_ecef(1);
                data.y_mat(i,j) = pos_ecef(2);
                data.z_mat(i,j) = pos_ecef(3);
            end
        end
    end
else % Grid was passed as input argument
    data.x_mat = grid_params.grid(:,:,1);
    data.y_mat = grid_params.grid(:,:,2);
    data.z_mat = grid_params.grid(:,:,3);
end
if nargout>2, grid = cat(3, data.x_mat, data.y_mat, data.z_mat); end; % Return grid if requested

%% Call the backprojection function with the appropriate inputs
% Convert metadata from CPHD structure into structure that bpBasic function wants
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
% The following is the only line of actual IFP in this whole function:
complex_image = bpBasic(data);
% "Basebanding".  This centers the frequency support, even if the image is
% formed to a DEM.  When working on image domain data, we need an average
% center frequency across all pulses.
f_c          = mean(nbdata.SC0  + (nbdata.SCSS * ...
    double(ifp_params.sample_range(ceil(length(ifp_params.sample_range)/2))-1)));
c            = SPEED_OF_LIGHT; % Speed of light (m/s)
complex_image = complex_image.*...
    exp(-1i*4*pi*f_c/c*data.x_mat*grid_params.range_vect_hat(1)).*...
    exp(-1i*4*pi*f_c/c*data.y_mat*grid_params.range_vect_hat(2)).*...
    exp(-1i*4*pi*f_c/c*data.z_mat*grid_params.range_vect_hat(3));

%% Create SICD-like metadata structure.  A true SICD metadata structure is
% not really possible, since back-projection is not supported by SICD, but
% we will fill it in where possible.
if nargout>1
    ifp_params = setstructfields(ifp_params, grid_params);
    sicdmeta = bp_sicd_meta(cphd_meta,nbdata,ifp_params);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////