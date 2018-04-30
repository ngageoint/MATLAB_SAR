function [ output_data ] = deskewfile( input_filename, varargin )
%DESKEWFILE Performs deskew (centering of the spectrum on zero frequency) on a file of complex data
%    output_data = deskewfile(input_filename, 'PropertyName', PropertyValue, ...)
%
%       Property name     Description
% 
%       output_filename   Output data is stored to a complex single-
%                            precision float SIO file (with no metadata).
%                            Default is to not save to file.
%       block_size        Determines how many lines are processed in memory
%                            at once.  Set to Inf to process whole image at
%                            once.  Default is 500 lines.
%       dim               Dimension over which to perform deskew.  Default
%                            = 1 (slow-time).
%       azlimits          Min and max samples in azimuth over which to
%                            compute.  Default = full image.
%       rnglimits         Min and max lines in range over which to compute.
%                            Default = full image.
%       framenumber       Dataset to process if INPUT_FILENAME is a
%                            multi-image file.  Default is 1.
%
% If an output argument is given to the function call, the deskewed complex
% data will be returned in an in-memory MATLAB array to OUTPUT_DATA.
% WARNING: Make sure your computer has enough memory to hold resulting
% array, if you use this option.
%
% Implemented by: Wade Schwartzkopf, NGA/IDT
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
p.FunctionName = 'deskewfile';
p.parse(varargin{:});

save_to_sio=~any(strcmp(p.UsingDefaults,'output_filename'))||~isempty(p.Results.output_filename);
save_to_variable=(nargout>0);
if ~save_to_sio&&~save_to_variable % Nothing to do
    warning('DESKEWFILE:NO_OUTPUT','Neither output filename nor output variable was specified.  No ouput will be calculated.');
    return;
end

%% Open data file
reader_obj=open_reader(input_filename);
if iscell(reader_obj) % If file contains more than one image
    % Just use first image for now.  Later we can allow user to specify
    % specific images from multi-image files.
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
end
if any(strcmp(p.UsingDefaults,'rnglimits'))
    rnglimits=[1 meta.ImageData.NumRows];
else
    rnglimits=p.Results.rnglimits;
end

%% Calculate input parameters for deskewmem
[ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( meta, p.Results.dim );
if all(~DeltaKCOAPoly(:))
    warning('DESKEWFILE:NO_DELTAKCOAPOLY','DELTAKCOAPOLY is not defined or is zero.  No deskew being performed.');
    if nargout>0
        output_data = single(reader_obj.read_chip(azlimits,rnglimits));
    end
    reader_obj.close();
    return;
end

%% Setup loop for processing through image in blocks of azimuth lines (block_size=Inf will process entire image in memory)
if save_to_variable, output_data=zeros(diff(azlimits)+1,diff(rnglimits)+1); end;
if save_to_sio, writer_object=open_sio_writer(p.Results.output_filename); end;
first_row_in_block=rnglimits(1);
wb_hand=waitbar(0,'Deskewing data...');
while(first_row_in_block<rnglimits(2))
    row_indices=first_row_in_block:min(rnglimits(2),first_row_in_block+p.Results.block_size-1);
    deskewed_block = deskewmem(single(reader_obj.read_chip(azlimits,row_indices([1 end]))),...
        DeltaKCOAPoly,az_coords_m(azlimits(1):azlimits(2)),...
        rg_coords_m(row_indices),p.Results.dim,fft_sgn);
    if save_to_variable, output_data(:,row_indices-rnglimits(1)+1)=deskewed_block; end;
    if save_to_sio, writer_object.write_row(deskewed_block); end;
    first_row_in_block=row_indices(end)+1;
    waitbar(double(first_row_in_block-rnglimits(1))/double(diff(rnglimits)+1),wb_hand);
end

%% Cleanup
close(wb_hand);
reader_obj.close();
if save_to_sio, writer_object.close(); end;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////