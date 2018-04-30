function [ output_meta ] = meta2sicd_s1noise( domnode, meta_product )
%META2SICD_S1NOISE Converts Sentinel-1 noise XML file into SICD format
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Setup
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

% Compute the first line of each burst in overall TIFF file
first_line = zeros(numel(meta_product),1);
for i = 1:(numel(meta_product)-1)
    first_line(i+1) = first_line(i) + meta_product{i}.ImageData.NumCols;
end
% Read noise parameters from XML
num_noise_vecs=str2double(xp.evaluate(...
    'count(noise/noiseVectorList/noiseVector)',domnode));
noisevals = cell(num_noise_vecs,1);
pixel = cell(num_noise_vecs,1);
noise_line = zeros(num_noise_vecs,1);
for i = 1:num_noise_vecs
    noise_line(i) = str2double(xp.evaluate(...
        ['noise/noiseVectorList/noiseVector[' num2str(i) ']/line'],...
        domnode));
    pixel{i} = str2num(xp.evaluate(...
        ['noise/noiseVectorList/noiseVector[' num2str(i) ']/pixel'],...
        domnode));
    noisevals{i} = str2num(xp.evaluate(...
        ['noise/noiseVectorList/noiseVector[' num2str(i) ']/noiseLut'],...
        domnode));
end
mode = char(xp.evaluate('noise/adsHeader/mode',domnode));
fit_tolerance = 0.1; % Simple 1D data for IW should fit quite tightly
noisepoly_order = 6; % For IW, this fits quite tightly
if upper(mode(1))=='S'
    % Average all of the noise values across azimuth for a single 1D
    % polynomial that varies in range.  They are roughly similar.  A
    % precise solution would include a 2D fit.
    pixel = {cell2mat(pixel)};
    noisevals = {cell2mat(noisevals)};
    noise_line = 0;
    % Allow for more slop for SM since we have 2D noise values, but aren't
    % doing a full 2D fit.  Because we average across all azimuth, there is
    % not reason to attempt a tight fit, since it won't match across all
    % lines for any fit.
    fit_tolerance = 0.5;
    noisepoly_order = 3;
% Double check our assumptions about the noise data
elseif any(mod(noise_line(1:(end-1)), ...
        double(meta_product{1}.ImageData.NumCols))~=0)
    warning('META2SICD_S1NOISE:UNEXPECTED_LINE','Expected noise values at beginning of each burst.');
    return;
end
% Compute SICD noise polynomials
output_meta = cell(numel(meta_product),1);
for i = 1:numel(meta_product) % First and last are usually outside bursts
    % Which noise values should we use for this burst?  We use the ones
    % that correspond to the first column of the burst.
    line_ind = find(noise_line==first_line(i),1);
    coords_rg_m = (double(0:(meta_product{i}.ImageData.NumRows-1)) - ...
        double(meta_product{i}.ImageData.SCPPixel.Row)) * ...
        meta_product{i}.Grid.Row.SS;
    coords_rg_m = coords_rg_m(pixel{line_ind}(:)+1);
    old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
    noisepoly = polyfit(coords_rg_m(:), 10*log10(noisevals{line_ind}(:)), noisepoly_order);
    warning(old_state);
    if any(polyval(noisepoly, coords_rg_m(:)) - ...
            10*log10(noisevals{line_ind}(:)) > fit_tolerance)
        warning('META2SICD_S1NOISE:NOISE_TOLERANCE', ...
            'Noise model does not appear to fit data tightly.');    
    end
    % These are typically very close across all bursts, but we compute for
    % each bursts independently anyway.
    output_meta{i}.Radiometric.NoiseLevel.NoiseLevelType = 'ABSOLUTE';
    output_meta{i}.Radiometric.NoiseLevel.NoisePoly = noisepoly(end:-1:1).';
end
% Throw warning for old, potentially inaccurate data
last_az_time = xp.evaluate(...
        ['noise/noiseVectorList/noiseVector[' num2str(num_noise_vecs) ']/azimuthTime'],...
        domnode);
if datenum(char(last_az_time),'yyyy-mm-ddTHH:MM:SS')<datenum(2015, 11, 25)
    warning('META2SICD_S1NOISE:EARLY_DATE','Sentinel-1 radiometric data prior to November 25, 2015 might not be accurate.');    
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////