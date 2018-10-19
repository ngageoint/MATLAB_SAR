function convert_complex_data(inFilename, outFilename, varargin)
%CONVERT_COMPLEX_DATA Converts the specified input file (inFilename) to the
%specified output file (outFilename).  Additional conversion parameters
%are allowed as described below.
%
% Inputs:
%       inFilename:  Name of file to convert.  It must be of any type
%                    supported by the toolbox open_reader.
%       outFilename: The name of the output file
%
% Outputs:
%       None
%
% Allowed properties:
%       Property name         Description
%       -----------------------------------------------------------------
%       output_format         The output file format to write.  Allowed
%                             values are SICD, NITF, NRL, and SIO.  Default
%                             is SICD.
%       detect                Whether the output image should be remapped
%                             to an 8-bit detected (non-complex) image
%       frames                Set of each frame to convert.  Default is
%                             all.
%       showWaitbar           Whether a GUI wait bar will be
%                             displayed while the image conversion is
%                             progressing
%
% Written by: Tom Krauss and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open the input file and output files.  Note the input file must be of a
% format recognized by the toolbox open_reader.  Also, we'll only read the
% first image of multi-image files.
tempreader   = open_reader(inFilename);
if iscell(tempreader)
    reader = tempreader;
else
    reader = {tempreader};
end

p = inputParser;
p.addParamValue('output_format', 'sicd', @(x) any(strcmpi(x,{'sicd', 'nitf', 'sio', 'nrl'})))
p.addParamValue('pixeltype', '');
p.addParamValue('detect', false);
p.addParamValue('frames', 1:length(reader));
p.addParamValue('showWaitbar', true);
p.parse(varargin{:});

for i=1:length(p.Results.frames)
    if length(p.Results.frames)>1
        [path, name, ext] = fileparts(outFilename);
        number_format = '%.3d'; %.' num2str(log10(max(p.Results.frames))+1) 'd'];
        newOutFilename = fullfile(path,[name num2str(p.Results.frames(i),number_format) ext]);
    else
        newOutFilename = outFilename;
    end
    sicdmeta = reader{p.Results.frames(i)}.get_meta();
    
    % Pick the correct output image writer object
    if p.Results.showWaitbar, h=waitbar(0); end % Intialize waitbar
    if p.Results.detect
        % Detected, but still a slant-plane image with no ground projection
        constructor_args={'data_type', 'uint8', 'is_complex', false};
        if p.Results.showWaitbar, waitbar(0,h,'Computing image statistics...'); end
        datamean = estimatemean(reader{p.Results.frames(i)});
    elseif strcmp(p.Results.pixeltype,'RE16I_IM16I') && ...
            ~strcmp(sicdmeta.ImageData.PixelType,'RE16I_IM16I')
        % Convert float to int by scaling to full signed int16 range
        % Must update sicdmeta PixelType
        % It's ok to change this temp copy since its use below pertains to
        % both the input and output image ignoring PixelType.
        constructor_args={'data_type', 'int16'};
        sicdmeta.ImageData.PixelType = 'RE16I_IM16I';
        if p.Results.showWaitbar, waitbar(0,h,'Computing image statistics...'); end
        scl = scale_to_int16(reader{p.Results.frames(i)});
    else
        constructor_args={}; % Nothing to override.  Use values in SICD metadata structure
    end
    if strcmp(p.Results.output_format,'sicd')
        writer = SICDWriter(newOutFilename, sicdmeta);
    elseif strcmp(p.Results.output_format,'nitf')
        writer = NITFWriter(newOutFilename, sicdmeta, constructor_args{:});
    elseif strcmp(p.Results.output_format,'sio')
        writer = SIOWriter(newOutFilename, sicdmeta, constructor_args{:});
    elseif strcmp(p.Results.output_format,'nrl')
        writer = NRLWriter(newOutFilename, sicdmeta, constructor_args{:});
    end
    
    
    % We'll write the data chunk-by-chunk.  Each chunk is a set of
    % consecutive lines.
    MAX_BLOCK_SIZE = 2^26; % in bytes (roughly 50MB)
    rowsPerBlock=max(1,floor(MAX_BLOCK_SIZE/(double(sicdmeta.ImageData.NumCols)*8))); % Assume 8-byte elements (largest allowed in SICD)
    numBlocks = ceil(double(sicdmeta.ImageData.NumRows)/rowsPerBlock);
    for j = 1:numBlocks
        blockStart = 1 + (j-1)*rowsPerBlock;
        blockEnd   = min( blockStart+rowsPerBlock-1, sicdmeta.ImageData.NumRows );
        if p.Results.showWaitbar
            waitbar(double(blockStart)/double(sicdmeta.ImageData.NumRows),h,'Writing data...');
        end
        
        data = reader{p.Results.frames(i)}.read_chip([1 sicdmeta.ImageData.NumCols],[blockStart blockEnd]);
        if p.Results.detect
            writer.write_chip(amplitudetodensity(data,30,40,datamean), [1 blockStart]);
        elseif exist('scl','var')
            writer.write_chip(int16(data * scl), [1 blockStart]);
        else
            writer.write_chip(data, [1 blockStart]);
        end
    end
    if p.Results.showWaitbar, close(h); end;
end

reader{1}.close(); % Should only need to close one, since only one file
% Although not necessary to call explicitly, the writer object will do
% footer writing, file closing, and other cleanup upon leaving scope.
% writer.delete();

end


function sample_mean = estimatemean(reader)
    sicdmeta = reader.get_meta();
    samplesize=[1000 1000]; % Exract 1000x1000 array of samples to estimate mean
    subsample=ceil(double([sicdmeta.ImageData.NumCols sicdmeta.ImageData.NumRows])./samplesize);
    data = abs(single(reader.read_chip([1 sicdmeta.ImageData.NumCols],...
            [1 sicdmeta.ImageData.NumRows],subsample)));
    sample_mean=mean(data(isfinite(data(:))));
end


% Use full int16 range and do not clip!
% (else autofocus, apodization etc. can break)
% Should this subsample like sample_mean() to avoid requiring full image in memory?
function scl = scale_to_int16(reader)
    sicdmeta = reader.get_meta();
    data = reader.read_chip([1 sicdmeta.ImageData.NumCols],...
                            [1 sicdmeta.ImageData.NumRows]);
    ifd = isfinite(data(:));
    scl = 32767 / max(max(abs(real(data(ifd)))), max(abs(imag(data(ifd)))));
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////