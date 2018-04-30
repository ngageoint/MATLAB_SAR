function [ readerobj ] = open_tiff_reader(filename)
%OPEN_TIFF_READER Intiates a reader object for TIFF imagery
%
% Works for most generic TIFFs, but really intended for SAR datasets stored
% in TIFF format.  For some GeoTIFFs (like RADARSAT-2 and Sentinel-1) this
% will also find and process metadata into the SICD format.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Sensor-specific metadata
% A TIFF file by itself is usually light on any SAR metadata.  Here we look
% for sensor-specific metadata files that might be associated with this
% TIFF file and use them to compute a SICD metadata structure.
tiff_sensor = guess_tiff_sensor(filename);
meta2sicd_tiff_fun = ['meta2sicd_tiff' tiff_sensor];
if ~isempty(tiff_sensor) && exist(meta2sicd_tiff_fun,'file')==2 % Sensor-specific
    [meta, symmetry] = feval(meta2sicd_tiff_fun, filename);
else
    meta = struct();
    symmetry = [0 0 0]; % Default orientation if none can be determined
end
    
%% Open reader
tiff_meta = imfinfo(filename);
if isfield(tiff_meta,'SampleFormat') && ... % Normal MATLAB TIFF reader does not handle complex data types
        any(strncmpi(tiff_meta.SampleFormat,'complex',7))
    symmetry(3) = ~symmetry(3); % Our TIFF reader is transposed from MATLAB's
    if iscell(meta) % TIFF should be split into multiple parts
        % This case really only works for Sentinel1.  We should make a more
        % generic solution that allows the meta2sicd function to define how
        % TIFF is split.
        last_col = 0;
        base_reader = open_ctiff_reader_noxml(filename,symmetry);
        for i = 1:numel(meta)
            readerobj{i} = subset_reader(base_reader,...
                    last_col + [1 meta{i}.ImageData.NumCols],[1 meta{i}.ImageData.NumRows]);
            last_col = last_col + meta{i}.ImageData.NumCols;
        end
    else
        readerobj = open_ctiff_reader_noxml(filename,symmetry);
    end
else
    readerobj = open_tiff_reader_noxml(filename,symmetry);
end
if iscell(meta)
    for i = 1:numel(meta)
        meta{i} = setstructfields(meta{i},readerobj{i}.get_meta()); % Merge metadata
        readerobj{i}.get_meta=@() meta{i};
    end
else
    meta = setstructfields(meta,readerobj.get_meta()); % Merge metadata
    readerobj.get_meta=@() meta;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////