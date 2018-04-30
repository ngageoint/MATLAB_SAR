function [ output_data ] = normalize_complex_file( input_filename, varargin )
%NORMALIZE_COMPLEX_FILE Normalize complex SAR data to make it have consistent
%characteristics in frequency domain
%
% Normalized complex data will have
%    1) Centered and deskewed frequency support
%    2) Uniform weighting
%    3) FFT sign of -1
%
% This function only applies these transforms in one dimension since we
% normally only process one dimension at a time, and since in general,
% frequency support deskew can only be applied in one dimension.
%
%    output_data = normalize_complex_file(input_filename, 'PropertyName', PropertyValue, ...)
%
%       Property name     Description
% 
%       output_filename   Output data is stored to a complex single-
%                            precision float SIO file.  Default is to not
%                            save to file.
%       block_size        Determines how many lines are processed in memory
%                            at once.  Set to Inf to process whole image at
%                            once.  Default is 500 lines.
%       dim               Dimension over which to perform normalization.
%                            Default = 1 (slow-time).
%       azlimits          Min and max samples in azimuth over which to
%                            compute.  Default = full image.
%       rnglimits         Min and max lines in range over which to compute.
%                            Default = full image.
%       framenumber       Dataset to process if INPUT_FILENAME is a
%                            multi-image file.  Default is 1.
%
% If an output argument is given to the function call, the normalized
% complex data will be returned in an in-memory MATLAB array to
% OUTPUT_DATA. WARNING: Make sure your computer has enough memory to hold
% resulting array, if you use this option.
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('output_filename','');
p.addParamValue('block_size',500);
p.addParamValue('dim',1);
p.addParamValue('azlimits',[]);
p.addParamValue('rnglimits',[]);
p.addParamValue('framenumber',1);
p.FunctionName = 'normalize_complex_file';
p.parse(varargin{:});

save_to_sio=~any(strcmp(p.UsingDefaults,'output_filename'))||~isempty(p.Results.output_filename);
save_to_variable=(nargout>0);
if ~save_to_sio&&~save_to_variable % Nothing to do
    warning('NORMALIZE_COMPLEX_FILE:NO_OUTPUT','Neither output filename nor output variable was specified.  No ouput will be calculated.');
    return;
end

%% Open data file
reader_obj=open_reader(input_filename);
if iscell(reader_obj) % If file contains more than one image
    for i=1:length(reader_obj)
        if i~=p.Results.framenumber
            reader_obj{i}.close();
        end
    end
    reader_obj=reader_obj{p.Results.framenumber};
end
meta=reader_obj.get_meta();
% Calculate AOI
if any(strcmp(p.UsingDefaults,'azlimits'))
    azlimits=[1 meta.ImageData.NumCols];
else
    azlimits=p.Results.azlimits;
    if strcmp(azlimits,'full')
        azlimits=[1 meta.ImageData.NumCols];
    end
end
if any(strcmp(p.UsingDefaults,'rnglimits'))
    rnglimits=[1 meta.ImageData.NumRows];
else
    rnglimits=p.Results.rnglimits;
    if strcmp(rnglimits,'full')
        rnglimits=[1 meta.ImageData.NumRows];
    end
end

%% Calculate input parameters for normalize_complex_mem
if p.Results.dim==1
    sicd_grid_struct = meta.Grid.Col;
else
    sicd_grid_struct = meta.Grid.Row;
end
% Vectors describing range and azimuth distances from SCP (in meters) for rows and columns
[ DeltaKCOAPoly, az_coords_m, rg_coords_m ] = deskewparams( meta, p.Results.dim );
coords_m_args = {az_coords_m, rg_coords_m};
coords_ind = {azlimits(1):azlimits(end), rnglimits(1):rnglimits(end)};
other_dim = 3-p.Results.dim;
coords_m = coords_m_args{other_dim};
coords_m_args{p.Results.dim}=coords_m_args{p.Results.dim}(coords_ind{p.Results.dim});
% Zeropad
if isfield(sicd_grid_struct,'ImpRespBW')&&isfield(sicd_grid_struct,'SS')
    zeropad = max(1,1/(sicd_grid_struct.ImpRespBW*sicd_grid_struct.SS));
else
    zeropad = 1;
end
% Weighting function
try
    weight_fun = sicdweight2fun(sicd_grid_struct);
catch % Unable to determine weighting function from metadata.
    % Estimate weighting from complex data instead
    weight_fun = estimate_weighting_file(reader_obj, p.Results.dim);
end
% FFT sign
if isfield(sicd_grid_struct,'Sgn')
    fft_sign = sicd_grid_struct.Sgn;
else
    fft_sign = -1;
end

%% Setup loop for processing through image in blocks of azimuth lines
if save_to_variable
    output_data = zeros(diff(azlimits)+1,diff(rnglimits)+1);
    output_data_args = {1:size(output_data,1),1:size(output_data,2)};
end
if save_to_sio
    writer_object=open_sio_writer(p.Results.output_filename);
    if p.Results.dim==1
        writer_object.writeline=writer_object.write_row;
    else
        writer_object.writeline=writer_object.write_column;
    end
end
if p.Results.dim==1
    linelimits = rnglimits;
    read_lines_fun=@(lims) reader_obj.read_chip(azlimits,lims);
else
    linelimits = azlimits;
    read_lines_fun=@(lims) reader_obj.read_chip(lims,rnglimits);
end
wb_hand=waitbar(0,'Normalizing data...');
first_line_in_block = linelimits(1);
while(first_line_in_block<linelimits(2))
    line_indices=first_line_in_block:...
        min(linelimits(2), first_line_in_block + p.Results.block_size - 1);
    coords_m_args{other_dim} = coords_m(line_indices);
    normalized_block = normalize_complex_mem(...
        single(read_lines_fun(line_indices([1 end]))), DeltaKCOAPoly, ...
        coords_m_args{:}, weight_fun, zeropad, fft_sign, p.Results.dim);
    if save_to_variable
        output_data_args{other_dim} = line_indices - linelimits(1) + 1;
        output_data(output_data_args{:}) = normalized_block;
    end
    if save_to_sio
        writer_object.writeline(normalized_block);
    end
    first_line_in_block = line_indices(end) + 1;
    waitbar(double(first_line_in_block-linelimits(1))/double(diff(linelimits)+1),wb_hand);
end

%% Cleanup
close(wb_hand);
reader_obj.close();
if save_to_sio, writer_object.close(); end;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////