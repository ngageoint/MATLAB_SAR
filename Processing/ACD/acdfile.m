function full_acd = acdfile(referenceImageName, matchImageName, startPoint_1, chipSize, skip, varargin)
%ACDFILE Creates a colorized amplitude change detection image from two overlapping images.
%
% Match image is interpolated to ground lat/long positions of the pixel
% from the reference image.  The reference image is then put into the red
% channel while the match image is put into "cyan".
%
% Inputs:
%       referenceImageName: Name of the reference image
%       matchImageName:     Name of the match image
%       startPoint_1:       The starting pixel location in the reference
%                           image (two-vector of row,column)
%       chipSize:           The size of the extracted chip (two-vector of
%                           width, height)
%       skip:               The "stride" of pixels to extract within the
%                           specified region.  A skip of 1 will take every
%                           pixel, a skip of 2 will take every other pixel,
%                           etc.
%       projectionType:     Either 'slant' (registers match image to
%                           reference slant plane) or 'ground_northup'
%                           (registers both images to a regular lat/lon
%                           grid).  Default is 'slant'.
%
% Outputs:
%       full_acd:           A full color image of the specified region with
%                           the specified stride (skip).
%
% Warning: AOI selected must be small enough to fit into memory.
%
% Written by: Tom Krauss, NGA/IDT
% Modified by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Parse input parameters
p = inputParser;
p.KeepUnmatched = true;
 % Default is to leave in slant plane of reference image
p.addParamValue('projectionType', 'slant', @(x) any(strcmp(x,{'slant','ground_northup'})));
p.addParamValue('applySVA', false);
p.parse(varargin{:});

% Get the metadata for both images
fprintf(1,'\n');
fprintf(1,'Reading imagery...');
tic;
reader_1     = open_reader(referenceImageName);
meta_1       = reader_1.get_meta();

reader_2     = open_reader(matchImageName);
meta_2       = reader_2.get_meta();

% Make sure that chipped section from match image contains all 4
% corners (and thus the entire AOI) from reference image.
corners=zeros(2,4);
corners(:,1) = startPoint_1;
corners(:,2) = startPoint_1+[chipSize(1)-1 0];
corners(:,3) = startPoint_1+[0 chipSize(2)-1];
corners(:,4) = startPoint_1+chipSize-1;
endPoint_1=corners(:,4);
pos_lla       = point_slant_to_ground(flipud(corners), meta_1);
corners_match = flipud(point_ground_to_slant(pos_lla, meta_2));
startPoint_2=min(corners_match,[],2);
endPoint_2=max(corners_match,[],2);
% Make sure match AOI indices are within image bounds
startPoint_2=min(double([meta_2.ImageData.NumCols; meta_2.ImageData.NumRows]),...
    max(1,startPoint_2));
endPoint_2=min(double([meta_2.ImageData.NumCols; meta_2.ImageData.NumRows]),...
    max(1,endPoint_2));

% Read AOI data from image files
oldImage_1 = double(reader_1.read_chip([startPoint_1(1) endPoint_1(1)], ...
    [startPoint_1(2) endPoint_1(2)], ...
    skip));
oldImage_2 = double(reader_2.read_chip([startPoint_2(1) endPoint_2(1)], ...
    [startPoint_2(2) endPoint_2(2)], ...
    skip));
T = toc;
fprintf('%.1f seconds\n', T)

% Upsample and apply SVA
if (p.Results.applySVA)
    % Upsample the image to the next integer value of "oversample factor".
    row_oversample = 1/(meta_1.Grid.Row.SS*meta_1.Grid.Row.ImpRespBW);
    col_oversample = 1/(meta_1.Grid.Col.SS*meta_1.Grid.Col.ImpRespBW);
    new_row_rate   = ceil(row_oversample);
    new_col_rate   = ceil(col_oversample);
    oldImage_1 = UpsampleImage(oldImage_1, ...
                               new_row_rate/row_oversample, ...
                               new_col_rate/col_oversample);
    oldImage_1 = Apodize2D_IIQ(oldImage_1, new_row_rate, new_col_rate);
    oldImage_1 = UpsampleImage(oldImage_1, ...
                               row_oversample/new_row_rate, ...
                               col_oversample/new_col_rate);

    row_oversample = 1/(meta_2.Grid.Row.SS*meta_2.Grid.Row.ImpRespBW);
    col_oversample = 1/(meta_2.Grid.Col.SS*meta_2.Grid.Col.ImpRespBW);
    new_row_rate   = ceil(row_oversample);
    new_col_rate   = ceil(col_oversample);
    oldImage_2 = UpsampleImage(oldImage_2, ...
                               new_row_rate/row_oversample, ...
                               new_col_rate/col_oversample);
    oldImage_2 = Apodize2D_IIQ(oldImage_2, new_row_rate, new_col_rate);
    oldImage_2 = UpsampleImage(oldImage_2, ...
                               row_oversample/new_row_rate, ...
                               col_oversample/new_col_rate);
end

% Calculate ground positions of each pixel in the two images
fprintf('Projecting images...');
tic;
[lats_1, lons_1] = project_image(startPoint_1, size(oldImage_1), skip, meta_1, varargin{:});
[lats_2, lons_2] = project_image(startPoint_2, size(oldImage_2), skip, meta_2, varargin{:});
T = toc;
fprintf('%.1f seconds\n', T)


% Create grid to which to interpolate
if strcmp(p.Results.projectionType,'ground_northup')
    % Uniform lat/lon grid to interpolate to.  Here we're using the same
    % number of output samples (it should probably be slightly bigger since
    % the image area is probably larger).
    pos_lla2     = point_slant_to_ground(flipud(corners_match), meta_2);
    chip_rows = size(oldImage_1,1);%%% * (max(pos_lla2(1,:))-min(pos_lla2(1,:)))/(max(pos_lla(1,:)) - min(pos_lla(1,:)));
    chip_cols = size(oldImage_1,2);%%% * (max(pos_lla2(2,:))-min(pos_lla2(2,:)))/(max(pos_lla(2,:)) - min(pos_lla(2,:)));
    [qx,qy] = meshgrid(linspace(min(pos_lla2(1,:)), max(pos_lla2(1,:)), chip_rows),...
                       linspace(min(pos_lla2(2,:)), max(pos_lla2(2,:)), chip_cols));
else % Remain in slant plane of reference image
    qx=lats_1; qy=lons_1;
end

% Build "Scattered data interpolant".   Delaunay triangulation of
% non-uniform gridding
fprintf(1, 'Interpolating images...');
tic;
% F_2 = TriScatteredInterp(lats_2(:), lons_2(:), abs(oldImage_2(:)) );
% outputImage_2 = F_2(qx,qy);
F_2_r = TriScatteredInterp(lats_2(:), lons_2(:), real(oldImage_2(:)) );
F_2_i = TriScatteredInterp(lats_2(:), lons_2(:), imag(oldImage_2(:)) );

re = F_2_r(qx,qy);
idx = isnan(re);
re(idx) = 0.0;
im = F_2_i(qx,qy);
idx = isnan(im);
im(idx) = 0.0;
outputImage_2 = abs(re + 1i*im);
if strcmp(p.Results.projectionType,'ground_northup') % First image must also be interpolated
    % F_1 = TriScatteredInterp(lats_1(:), lons_1(:), abs(oldImage_1(:)) );
    % outputImage_1 = F_1(qx,qy);
    F_1_r = TriScatteredInterp(lats_1(:), lons_1(:), real(oldImage_1(:)) );
    F_1_i = TriScatteredInterp(lats_1(:), lons_1(:), imag(oldImage_1(:)) );
    re = F_1_r(qx,qy);
    idx = isnan(re);
    re(idx) = 0.0;
    im = F_1_i(qx,qy);
    idx = isnan(im);
    im(idx) = 0.0;
    outputImage_1 = abs(re + 1i*im);
    outputImage_1 = fliplr(outputImage_1); % fliplr orients north-up
    outputImage_2 = fliplr(outputImage_2); % image will be transposed at display
else
    outputImage_1 = abs(oldImage_1); % Reference image remains the same
end
T = toc;
fprintf(1, '%.1f seconds\n', T);


% For points outside of the image the interpolator (bilinear) will fail.
% The failure is indicated by a "NaN" value.  We'll find these failed
% regions (pixels) and just set them to zero making the image area black.
% Note, we're actually using a small non-zero number because the remapping
% done later takes a log.
outputImage_2(isnan(outputImage_2)) = eps;
if strcmp(p.Results.projectionType,'ground_northup')
    outputImage_1(isnan(outputImage_1)) = eps;
end

% After we have intepolated images, do a basic (translation only) sub-pixel
% FFT-based registration.
aoi_buffer=32; % Estimate of how far geolocation might be off
% DFT registration "wraps" around.  If no zeropadded, wraparound will
% overlap good image data.  Remove zeropad after dft registration.
outputImage_1=padarray(outputImage_1,[aoi_buffer aoi_buffer]);
outputImage_2=padarray(outputImage_2,[aoi_buffer aoi_buffer]);
[outputStats outputImageSpect] = dftregistration(fft2(outputImage_1),fft2(outputImage_2),20);
outputImage_2 = abs(ifft2(outputImageSpect));
outputImage_1=outputImage_1((aoi_buffer+1):(end-aoi_buffer),(aoi_buffer+1):(end-aoi_buffer));
outputImage_2=outputImage_2((aoi_buffer+1):(end-aoi_buffer),(aoi_buffer+1):(end-aoi_buffer));

% Combine the channels to make RGB red/cyan ACD image.
full_acd = acdmem(outputImage_1, outputImage_2);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////