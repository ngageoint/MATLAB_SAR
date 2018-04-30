function complex_image = bpBasic(data)
%BPBASIC This function performs a basic backprojection operation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following fields need to be populated:                           %
% data.deltaF:  Vector containing the step size of frequency data (Hz) %
% data.minF:  Vector containing the start frequency of each pulse (Hz) %
% data.x_mat:  The x-position of each pixel (m)                        %
% data.y_mat:  The y-position of each pixel (m)                        %
% data.z_mat:  The z-position of each pixel (m)                        %
% data.Tx.X:  The transmit x-position of the sensor at each pulse (m)  %
% data.Tx.Y:  The transmit y-position of the sensor at each pulse (m)  %
% data.Tx.Z:  The transmit z-position of the sensor at each pulse (m)  %
% data.R0:  The (bistatic) range to motion-comp point (m)              %
% data.phdata:  Phase history data (frequency domain)                  %
%               Fast time in rows, slow time in columns                %
%                                                                      %
% The following fields are optional:                                   %
% data.Rcv.X:  The receive x-position of the sensor at each pulse (m)  %
% data.Rcv.Y:  The receive y-position of the sensor at each pulse (m)  %
% data.Rcv.Z:  The receive z-position of the sensor at each pulse (m)  %
%              If no Rcv field is given, the monostataic case          %
%              (Rcv = Tx) is assumed.                                  %
% data.Nfft:  Size of the FFT to form the range profile                %
% data.hide_waitbar: Whether to suppress display of waitbar or not.    %
%                    If you call bpBasic from within a loop, you will  %
%                    probably want to use this.  Default is false.     %
%                                                                      %
% The output is:                                                       %
% complex_image:  The complex image value at each pixel                %
%                                                                      %
% History:                                                             %
%    Original code                                                     %
%       LeRoy Gorham, Air Force Research Laboratory, WPAFB, OH         %
%       leroy.gorham@wpafb.af.mil                                      %
%       Original Date Released:  8 Apr 2010                            %
%    Modified by Wade Schwartzkopf, NGA/IDT                            %
%       Wade.C.Schwartzkopf.ctr@nga.mil                                %
%       Structural changes to handle CPHD                              %
%    Modified by Daniel Andre, Dstl, Porton Down, UK                   %
%       dandre@dstl.gov.uk                                             %
%       Extension of the code to the bistatic case (1 Sep 2011)        %
%                                                                      %
% Gorham, L.A. and Moore, L.J., "SAR image formation toolbox for       %
%   MATLAB,"  Algorithms for Synthetic Aperture Radar Imagery XVII     %
%   7669, SPIE (2010).                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup input parameters
% Determine the size of the phase history data for easier to read code
num_pulses = size(data.phdata,2);

% Maximum scene extent in range (m)
c = 299792458; % Speed of light (m/s)
max_extent_range=c/mean(data.deltaF); % Two-way range (transmitter, target, receiver)

% Setup reasonable default value for Nfft if none passed in
if ~isfield(data,'Nfft')
    data.Nfft = 2^(3+nextpow2(size(data.phdata,1))); % Use power-of-2 FFT
    % The linear interolation in the interp1 below causes cross-range
    % artifacts unless we sinc interpolate first with a much larger FFT size
    % than necessary.  Thus the "3+" in the line above.
end

% Assume monostatic case if bistatic info is not given
if ~isfield(data,'Rcv')
    data.Rcv = data.Tx;
end

show_waitbar = ~(isfield(data,'hide_waitbar')&&data.hide_waitbar);

%% Do actual backprojection computation
% Calculate the range to every bin in the range profile (m)
range_vector = linspace(-data.Nfft/2,data.Nfft/2-1,data.Nfft)*max_extent_range/data.Nfft;

% Initialize the image with all zero values
complex_image = zeros(size(data.x_mat));

if show_waitbar
    wb = waitbar(0);
    tic;
end
% Loop through every pulse
for ii = 1:num_pulses
    % Form the range profile with zero padding added
    rc = fftshift(ifft(data.phdata(:,ii),data.Nfft));

    % Calculate differential (bistatic) range for each pixel in the image (m)
    dR = sqrt((data.Tx.X(ii)-data.x_mat).^2 + ... % Range from transmit to pixel
              (data.Tx.Y(ii)-data.y_mat).^2 + ...
              (data.Tx.Z(ii)-data.z_mat).^2) + ...
         sqrt((data.Rcv.X(ii)-data.x_mat).^2 + ... % Range from receive to pixel
              (data.Rcv.Y(ii)-data.y_mat).^2 + ...
              (data.Rcv.Z(ii)-data.z_mat).^2) - ...
         (2*data.R0(ii)); % Range from transmit to motion-comp point to receive

    % Calculate phase correction for image
    phCorr = exp(1i*2*pi*data.minF(ii)/c*dR);

    % Determine which pixels fall within the range swath
    I = find(and(dR > min(range_vector), dR < max(range_vector)));

    % Update the image using linear interpolation
    complex_image(I) = complex_image(I) + interp1(range_vector,rc,dR(I),'linear') .* phCorr(I);
            
    if show_waitbar
        % Determine remaining execution time and display
        t_sofar = toc;
        t_est = (t_sofar*num_pulses/ii)-t_sofar;
        wb_message=sprintf('Pulse %d of %d, Time remaining: %s',ii,...
            num_pulses,datestr(datenum(0,0,0,0,0,t_est),13));
        waitbar(ii/num_pulses,wb,wb_message);
    end
end
if show_waitbar
    close(wb);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////