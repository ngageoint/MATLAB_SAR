function [ xpoly, ypoly, zpoly, avg_resid ] = sv2poly( xpos, ypos, zpos, ...
    times, xvel, yvel, zvel )
%SV2POLY Compute best fit polynomials from 3D state vectors
%
% This function determines a polynomial order that achieves good precision
% without being unreasonably high or overfitting.  One would generally only
% use this function on state vectors from an orbital system, where the
% positions/velocity would reasonably follow a smooth path.
%
% "Optimal" polynomial order depends on both the length of the collect
% around its orbit and the noise in the state vector data.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Check arguments
if ~isequal(size(xpos), size(ypos), size(zpos), size(times)) || ...
        (exist('zvel','var') && ...
        ~isequal(size(xvel), size(yvel), size(zvel), size(times)))
    error('SV2POLY:INVALID_INPUT_ARGS','All state vector components must be same size/shape.');
end
% Vectorize inputs to assure column vectors
times = times(:);
xpos = xpos(:);
ypos = ypos(:);
zpos = zpos(:);
if exist('zvel', 'var')
    xvel = xvel(:);
    yvel = yvel(:);
    zvel = zvel(:);
end

% Because of the scale of the values used in orbital state vectors,
% MATLAB's polyfit oftens detects them as badly conditioned, but results
% work out OK anyway for reasonable polynomail orders. We will assure in
% this function that polynomial order is reasonable.
old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

%% Method 1.  Naive approach.
% Hardcode polynomial order.  This is usually sufficient.  5th order
% generally works well for most spaceborne SAR collects of typical length
% (< 200 s), often putting evaluation of the polynomials well within 1 mm
% of error.
% polyorder = min(5, numel(times) - 1);

if exist('zvel', 'var')
    %% Method 2.  Cross check with velocity.
    % One could also find the order of polynomial that most accurately
    % describes this position, but use velocity (which is presumably
    % measured independently) as cross-validation so that the data is not
    % being overfit.
    for polyorder = 1:(numel(times)-1)
        xpoly = polyfit(times, xpos, polyorder);
        ypoly = polyfit(times, ypos, polyorder);
        zpoly = polyfit(times, zpos, polyorder);
        error_list(polyorder) = mean(sqrt(sum(([xvel, yvel, zvel] - ...
            [polyval(polyder(xpoly), times), ...
            polyval(polyder(ypoly), times), ...
            polyval(polyder(zpoly), times)]).^2, 2)));
    end
    % Increase order only as long as it results in a "significant" decrease
    % (half) in the error of the velocity.
    polyorder = find(error_list(1:(end-1))<(2*error_list(2:end)), 1, 'first');
    if isempty(polyorder), polyorder = numel(error_list); end  % All orders improve quality
else
    %% Method 3. "Significant" improvement.
    % One could also use the lowest order of polynomial that significantly
    % reduces the error of the position fit.
    for polyorder = 1:(numel(times)-2)  % polyorder of n-1 will always be exact match
        % mu term used for error computation in loop since we are
        % potentially computing high order polynomials, generally well
        % beyond what is required for a good fit, and these will likely be
        % badly conditioned. For our final fit, after the polynomial order
        % is determined, we will assume order and fit is reasonable-- even
        % if MATLAB would have thrown a warning.
        [xpoly, ~, xmu] = polyfit(times, xpos, polyorder);
        [ypoly, ~, ymu] = polyfit(times, ypos, polyorder);
        [zpoly, ~, zmu] = polyfit(times, zpos, polyorder);
        error_list(polyorder) = mean(sqrt(sum(([xpos, ypos, zpos] - ...
            [polyval(xpoly, times, [], xmu), ...
            polyval(ypoly, times, [], ymu), ...
            polyval(zpoly, times, [], zmu)]).^2, 2)));
    end
    % Increasing polynomial order by one must result in error being cut in
    % half in order to be considered "significant".
    polyorder = find(error_list(1:(end-1))<(2*error_list(2:end)), 1, 'first');
    if isempty(polyorder), polyorder = numel(error_list); end  % All orders improve quality
end

% Once optimal polynomial order is determined, do actual polynomial fit
xpoly  = polyfit(times, xpos, polyorder);
ypoly  = polyfit(times, ypos, polyorder);
zpoly  = polyfit(times, zpos, polyorder);

warning(old_state);

% Compute final residuals if requested
if nargout>3
    avg_resid = mean(sqrt(sum(([xpos, ypos, zpos] - ...
        [polyval(xpoly, times), ...
        polyval(ypoly, times), ...
        polyval(zpoly, times)]).^2, 2)));
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////