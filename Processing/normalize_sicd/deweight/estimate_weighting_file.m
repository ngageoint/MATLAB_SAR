function [ weight_fun ] = estimate_weighting_file( filename, dim )
%ESTIMATE_WEIGHTING_FILE Estimate the weighting applied to complex SAR data file
%
%    weight_fun = estimate_weighting_file(filename, dimension)
%
%       Parameter name    Description
% 
%       filename          Filename of complex SAR data file.  Can also be a
%                         reader object from open_reader.m, rather than a
%                         filename.
%       dim               Dimension over which to estimate weighting.  1 =
%                         azimuth, 2 = range.  Default = 1;
%       weight_fun        Function handle that generates weighting.  Takes
%                            a single input parameter, which is the number
%                            of elements in the resulting weighting vector.
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('dim','var')
    dim = 1;
end
other_dim = 3-dim;

%% Metadata parsing
if ischar(filename) && exist(filename,'file') % filename is string
    reader_object = open_reader(filename);
    close_reader = true;
else % filename is open_reader.m reader object
    reader_object = filename;
    close_reader = false;
end
% This section is a lot of cut-and-paste redundant code from
% normalize_complex_file.m.  Surely there is a cleaner way to do this...
meta = reader_object.get_meta();
try
    [ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( meta, dim );
    coords_m_args = {az_coords_m, rg_coords_m};
    coords_m = coords_m_args{other_dim};
catch % Insufficient metadata
    DeltaKCOAPoly = 0;
end
% Zeropad
if isfield(meta,'Grid')&&all(isfield(meta.Grid,{'Col','Row'}))
    if dim==1
        sicd_grid_struct = meta.Grid.Col;
    else
        sicd_grid_struct = meta.Grid.Row;
    end
    if isfield(sicd_grid_struct,'ImpRespBW')&&isfield(sicd_grid_struct,'SS')
        zeropad = max(1,1/(sicd_grid_struct.ImpRespBW*sicd_grid_struct.SS));
    else
        zeropad = 1;
    end
else
    zeropad = 1;
end

%% Sample data
% We will sample about a million pixels in full rows or columns spread
% across this image.  This should be enough to get a good sample, but not
% so many that it takes very long to read.
pix_to_sample = 1000000;
if dim==1
    lines_to_read = pix_to_sample/meta.ImageData.NumCols;
    lines_in_file = meta.ImageData.NumRows;
else
    lines_to_read = pix_to_sample/meta.ImageData.NumRows;
    lines_in_file = meta.ImageData.NumCols;
end
dec_factor = [1 1];
dec_factor(other_dim) = max(floor(lines_in_file/lines_to_read),1); % Sample full lines across image
data = single(reader_object.read_chip([],[],dec_factor));
if close_reader
    reader_object.close();
end

%% Estimate weighting
if any(DeltaKCOAPoly(:)~=0) % Center frequency support if required
    coords_m_args{other_dim} = coords_m(1:dec_factor(other_dim):end);
    data = deskewmem(data, DeltaKCOAPoly, coords_m_args{:}, dim, fft_sgn);
end
% estimate_weighting_mem assumes frequency support already centered
weight_fun = estimate_weighting_mem(data,dim,zeropad);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////