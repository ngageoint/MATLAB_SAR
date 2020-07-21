function [ readerobj ] = open_cnitf_reader( filename, imageind )
%OPEN_NITF_READER Intiates a reader for some types of NITF files
%
% Disclaimer: NITF is a very flexible format and has many different
% varieties.  The reader only handles a few types of NITF, mainly complex
% slant-plane imagery from a few sensors.  Hopefully the framework is there
% to include any NITF image data though, as the need arises for different
% varieties.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Read metdata from each image and convert to SICD
native_metadata=read_nitf_meta(filename);
if nargin<2, imageind=1:native_metadata.filehdr.NUMI; end;
for i=1:length(imageind)
    % Are there any sensor specific adjustments to the metadata we must make?
    nitf_format = guess_nitf_format(native_metadata.imagesubhdr{i});

    % Symmetry is required to properly define SICD fields and pixel order
    % We always try to orient for azimuth first dimension, view-from-above
    % (must do transpose to view shadows-down in MATLAB).  NITF puts some
    % orientation in the CMETAA IF_RANGE_DATA field.  There seems to be an
    % implicit "view-from-above" assumption in NITF, at least in all the
    % examples we have seen, so only rotations and not flips are allowed.
    symmetry_fun = ['symmetry_' nitf_format];
    if ~isempty(nitf_format) && exist(symmetry_fun,'file')==2 % Sensor-specific
        symmetry{i}=feval(symmetry_fun, native_metadata.imagesubhdr{i});
    elseif isfield(native_metadata.imagesubhdr{i},'CMETAA')
        % CMETAA IF_RANGE_DATA field is a generic sensor-independent way to
        % describe data orientation/symmetry
        switch native_metadata.imagesubhdr{i}.CMETAA.IF_RANGE_DATA
            case 'ROW_INC'
                symmetry{i}=[0 0 0];
            case 'ROW_DEC'
                symmetry{i}=[1 1 0]; % Rotate 180 degrees
            case 'COL_INC'
                symmetry{i}=[0 1 1]; % Need to rotate 90 degrees clockwise
            case 'COL_DEC'
                symmetry{i}=[1 0 1]; % Need to rotate 90 degrees counter-clockwise
            otherwise
                symmetry{i}=[0 0 0]; % Default orientation if none is given
        end
    else
        symmetry{i}=[1 1 0]; % Default orientation if none is given
    end
    % Convert metadata in the image subheader to SICD
    meta{i}=meta2sicd_nitf(native_metadata.imagesubhdr{i}, symmetry{i});
    % Make sensor specific adjustments to SICD metadata might be required
    metadata_fun = ['meta2sicd_' nitf_format];
    if ~isempty(nitf_format) && exist(metadata_fun,'file')==2
        % Some data providers do not fill out enough metadata to
        % sufficiently describe their data so we must make metadata
        % assumptions based on sensor type, or they fill out the
        % "sensor-independent" TREs incorrectly so we must override the
        % normally computed metadata, or they have their own
        % sensor-specific TREs that we must handle.
        meta{i} = feval(metadata_fun, meta{i}, native_metadata.imagesubhdr{i}, native_metadata);
    end
    % Compute derived SICD fields
    meta{i}=derived_sicd_fields(meta{i});
    if ~isfield(meta{i}, 'ImageData')||~isfield(meta{i}.ImageData,'PixelType')
        % Default unless we know otherwise, since we cast to this anyway
        meta{i}.ImageData.PixelType = 'RE32F_IM32F';
    end
    % Pass on native metadata so that it could be browsed as an extension
    % to SICD metadata structure
    % TREs in file extended header should be distributed to all images.
    if isfield(native_metadata.filehdr,'XHD')&&isstruct(native_metadata.filehdr.XHD)
        XHD_fieldnames = fieldnames(native_metadata.filehdr.XHD);
        for j=1:length(XHD_fieldnames)
            native_metadata.imagesubhdr{i}.(XHD_fieldnames{j}) = ...
                native_metadata.filehdr.XHD.(XHD_fieldnames{j});
        end
    end
    meta{i}.native.nitf=native_metadata.imagesubhdr{i};
end

%% Compute parameters necessary for reading image data
data_offset=zeros(native_metadata.filehdr.NUMI,1); % Calculate file offsets to the start of each image segment
nextImgSubhdroffset=native_metadata.filehdr.HL;
for i = 1:native_metadata.filehdr.NUMI
    data_offset(i)=nextImgSubhdroffset+native_metadata.filehdr.LISH(i);
    if isfield(native_metadata.imagesubhdr{i},'IMDATOFF')
        data_offset(i)=data_offset(i)+native_metadata.imagesubhdr{i}.IMDATOFF;
    end
    nextImgSubhdroffset=nextImgSubhdroffset+native_metadata.filehdr.LISH(i)+native_metadata.filehdr.LI(i);
end

%% Build object
readerobj={};
for i=1:length(imageind)
    datasize=[native_metadata.imagesubhdr{i}.NCOLS native_metadata.imagesubhdr{i}.NROWS];
    datatype=nitf2matlab_datatype(native_metadata.imagesubhdr{i});
    if native_metadata.imagesubhdr{i}.NBANDS>1
        band_order=[];
        for j = 1:native_metadata.imagesubhdr{i}.NBANDS
            band_order = [band_order native_metadata.imagesubhdr{i}.(['ISUBCAT' num2str(j)])];
        end
        if native_metadata.imagesubhdr{i}.NBANDS==2 && isempty(band_order)
            band_order = 'IQ'; % Default to this if no ISUBCAT given
        end
    end
    if native_metadata.imagesubhdr{i}.PVTYPE(1)~='C'&&... % Not complex
        native_metadata.imagesubhdr{i}.NBANDS==1
        complextype=0;
    elseif native_metadata.imagesubhdr{i}.NBANDS>1&&...
            all(band_order(1:2:end)=='I')&&all(band_order(2:2:end)=='Q')
        complextype=1;
    elseif native_metadata.imagesubhdr{i}.NBANDS==2&&...
            any(strcmp(band_order,{'MP','PM'}))
        if isfield(native_metadata.imagesubhdr{i},'CMETAA')
            mag_map=native_metadata.imagesubhdr{i}.CMETAA.CMPLX_MAG_REMAP_TYPE;
            lin_scale=native_metadata.imagesubhdr{i}.CMETAA.CMPLX_LIN_SCALE;
            linlog_tp=native_metadata.imagesubhdr{i}.CMETAA.CMPLX_LINLOG_TP;
        else
            mag_map='LINM';
            lin_scale=1;
            linlog_tp=0;
        end
        complextype=@(x) pm_to_complex(x, band_order, mag_map, lin_scale, linlog_tp) ;
    elseif native_metadata.imagesubhdr{i}.PVTYPE(1)=='C'&&... % Lockheed-Martin C32 format
            ((strcmp(native_metadata.imagesubhdr{i}.IC,'C0')&&...
            any(strcmp(native_metadata.imagesubhdr{i}.COMRAT,{'C32','03.1'})))||...
            strcmp(datatype,'float16'))
        datatype='uint16';
        complextype=@c32_to_complex;
    elseif native_metadata.imagesubhdr{i}.PVTYPE(1)=='C'&&...
            native_metadata.imagesubhdr{i}.IC(1)=='N'
        complextype=1;
    else
        error('OPEN_CNITF_READER:UnsupportedFormat','Unsupported NITF format');
    end
    endian='b'; % NITF always big-endian
    
    new_reader = {};
    if strcmpi(native_metadata.imagesubhdr{i}.IC,'C8')
        % J2K compressed files must be handled specially.  This is an UGLY
        % hack!  We save the J2K portion of the NITF to a standalone J2K
        % file, because that's the only way MATLAB knows how to read J2K.
        int_filename = tempname(); % Intermediate file
        fid = fopen(filename,'r','b'); % NITF file
        fseek(fid,data_offset(imageind(i)),'bof');
        fid2 = fopen(int_filename,'w','b'); % Standalone J2K file
        % Process data in blocks since it may not all fit into memory
        blocksize = 2^26; % About 64 MBytes
        num_blocks = ceil(native_metadata.filehdr.LI(imageind(i))/blocksize);
        wb = waitbar(0,'Processing J2K data');
        for j = 1:num_blocks
            j2k_data = fread(fid,blocksize,'uint8');
            fwrite(fid2,j2k_data,'uint8');
            waitbar(j/num_blocks,wb);
        end
        close(wb);
        clear j2k_data; % Big and we don't need it anymore
        fclose(fid);
        fclose(fid2);
        % For most NITF data, the order of the dimensions described in the
        % datasize variable depends on the order the pixels are written in
        % the file.  For this case, it refers to the order of the
        % dimensions in the MATLAB J2K API.
        datasize_sym = datasize;
        if symmetry{i}(3), datasize_sym=datasize(end:-1:1); end
        meta{i}.native.j2k = imfinfo(int_filename);
        levels = meta{i}.native.j2k.WaveletDecompositionLevels;
        new_reader=chipfun2readerobj(@(varargin) ...
            j2k_chipper(int_filename, datasize, symmetry{i}, levels, varargin{:}), datasize_sym);
        new_reader.close = @() delete(int_filename); % We will delete the intermediate file when we are done
    elseif (native_metadata.imagesubhdr{i}.IMODE=='P'&&...
            native_metadata.imagesubhdr{i}.NBPR==1)||...
            (native_metadata.imagesubhdr{i}.IMODE=='B'&&...
            native_metadata.imagesubhdr{i}.NBPR==1&&...
            native_metadata.imagesubhdr{i}.NBANDS==1)
        new_reader=open_generic_reader(filename, datasize,...
            datatype, complextype, data_offset(imageind(i)), endian, ...
            symmetry{i},1,datasize);
    elseif native_metadata.imagesubhdr{i}.IMODE=='B' &&...
            native_metadata.imagesubhdr{i}.NBANDS==1
        blocksize=[native_metadata.imagesubhdr{i}.NPPBH native_metadata.imagesubhdr{i}.NPPBV];
        new_reader=open_generic_reader(filename, datasize,...
            datatype, complextype, data_offset(imageind(i)), endian,...
            symmetry{i}, native_metadata.imagesubhdr{i}.NBANDS, blocksize);
    elseif native_metadata.imagesubhdr{i}.IMODE=='P'
        blocksize=[native_metadata.imagesubhdr{i}.NPPBH native_metadata.imagesubhdr{i}.NPPBV];
        bands = native_metadata.imagesubhdr{i}.NBANDS;
        if complextype && native_metadata.imagesubhdr{i}.PVTYPE(1)~='C'
            bands = bands/2;
        end
        new_reader=open_generic_reader(filename, datasize,...
            datatype, complextype, data_offset(imageind(i)), endian,...
            symmetry{i}, bands, blocksize);
    elseif native_metadata.imagesubhdr{i}.IMODE=='S'||...
            (native_metadata.imagesubhdr{i}.IMODE=='B'&&...
            native_metadata.imagesubhdr{i}.NBPR==1) % Only multi-band IMODE=B, since single-band caught above
        blocksize=[native_metadata.imagesubhdr{i}.NPPBH native_metadata.imagesubhdr{i}.NPPBV];
        band_offsets = prod(blocksize) * ... % Pixel per block
            native_metadata.imagesubhdr{i}.NBPR*native_metadata.imagesubhdr{i}.NBPC * ... % Total number of blocks
            (native_metadata.imagesubhdr{i}.NBPP/8) * ... % Pixel size in bytes
            (0:(native_metadata.imagesubhdr{i}.NBANDS-1));
        for j=1:native_metadata.imagesubhdr{i}.NBANDS
            [chipper_function{j}, close_function{j}]=generic_chipper( filename, datasize,...
                datatype, 0, data_offset(imageind(i))+band_offsets(j), endian, symmetry{i}, 1, blocksize);
        end
        [chipper_all_bands,close_all_bands]=build_mb_chipper(chipper_function, close_function, complextype);
        if symmetry{i}(3), datasize=datasize(end:-1:1); end
        new_reader = chipfun2readerobj(chipper_all_bands, datasize);
        new_reader.close=close_all_bands;
    % If all else fails resort to MATLAB's NITF reader
    elseif exist('nitfread') % Do we have the image processing toolbox?
        new_reader = chipfun2readerobj(@(varargin) ml_nitf_chipper(filename, datasize, complextype, symmetry{i}, i, varargin{:}), datasize);
        new_reader.close=@() 1; % Nothing to do
    else
        error('OPEN_CNITF_READER:UNSUPPORTED_NITF_TYPE','Unsupported NITF type.');
    end
    new_reader.get_meta=@() meta{i};
    readerobj{end+1} = new_reader;
end

if numel(readerobj)==1
    readerobj = readerobj{1};
end

end

% Multi-band chipper.  Handles band sequential by doing multiple reads and
% concatenating.  Waitbar behavior odd, since you get a waitbar for each
% band.
function [concatenated_chipper, concatenated_close] = build_mb_chipper(chippers, close_functions, complextype)
    concatenated_chipper=@mb_chipper;
    concatenated_close=@mb_close;
   
    function data_out = mb_chipper(varargin)
        data_out=chippers{1}(varargin{:});
        for i=2:length(chippers)
            data_out(:,:,i)=chippers{i}(varargin{:});
        end
        if isequal(complextype,0) % Not complex
            % Do nothing
        elseif isequal(complextype,1) % Standard I/Q
            data_out=complex(data_out(:,:,1:2:end),data_out(:,:,2:2:end));
        elseif isa(complextype,'function_handle') % Custom complex types
            data_out=complextype(data_out);
        end
    end

    function mb_close()
        for i=1:length(close_functions)
            close_functions{i}();
        end
    end
end

% MATLAB NITF chipper
% Uses the MATLAB built-in NITFREAD function.  This is limited and only
% works for a few types of NITF.
function data_out = ml_nitf_chipper(filename, datasize, complextype, symmetry, im_index, varargin)
    newargs=reorient_chipper_args(symmetry, datasize, varargin{:});
    rows=[1 1 datasize(2)];
    cols=[1 1 datasize(1)];
    if nargin>5&&~isempty(newargs{1})
        rows(1)=newargs{1}(1);
        rows(3)=newargs{1}(2);
    end
    if nargin>6&&~isempty(newargs{2})
        cols(1)=newargs{2}(1);
        cols(3)=newargs{2}(2);
    end
    if nargin>7&&~isempty(newargs{3})
        rows(2)=newargs{3}(1);
        cols(2)=newargs{3}(2);
    end
    data_out = nitfread(filename,im_index,'PixelRegion',{cols,rows});
    data_out = permute(data_out,[2 1 3]); % MATLAB convention different than NITF
    if isequal(complextype,0) % Not complex
        % Do nothing
    elseif isequal(complextype,1) % Standard I/Q
        data_out=complex(data_out(:,:,1:2:end),data_out(:,:,2:2:end));
    elseif isa(complextype,'function_handle') % Custom complex types
        data_out=complextype(data_out);
    end
    data_out=reorient_chipper_data(symmetry, data_out);
end

% Function for chipping the J2K intermediate file
function out = j2k_chipper(filename, datasize, symmetry, max_level, varargin)
    newargs=reorient_chipper_args(symmetry, datasize, varargin{:});
    [dim1range,dim2range,subsample]=check_chipper_args(datasize, newargs{:});
    % The J2K API doesn't allow for regular decimation as we do with most
    % of our readers.  So we read the nearest reduction level and resize to
    % the equivalent resolution of the requested decimation.  This will not
    % be binary equivalent to raw deceimation (maybe this is good, since
    % raw decimation results in poor quality), but visually it should
    % produce an acceptable reduced resolution image quickly.
    reduction_level = min(max_level, floor(log2(min(subsample))));
    out = imread(filename, 'j2c', 'PixelRegion', ...
        {(floor(([dim2range(1) dim2range(2)]-1)/(2^reduction_level))+1),...
        (floor(([dim1range(1) dim1range(2)]-1)/(2^reduction_level))+1)}, ...
        'ReductionLevel', reduction_level);
    out = out.'; % Convert to row-major order
    if any(subsample>1)
        out = imresize(out, [numel(dim1range(1):subsample(1):dim1range(2)) ...
            numel(dim2range(1):subsample(2):dim2range(2))]);
    end
    out=reorient_chipper_data(symmetry, out);
end

% Convert phase/magnitude data into complex
function data_out = pm_to_complex(data_in, order, mag_map, lin_scale_factor, ll_tp)
    % We do not bother with the CMETAA.CMPLX_PHASE_QUANT_FLAGS here since,
    % for most data we have worked with, it has always been NS.
    if ~exist('order','var'), order='MP'; end;
    if ~exist('mag_map','var'), mag_map='LINM'; end; % Default is no remap
    if ~exist('lin_scale_factor','var'), lin_scale_factor = 1; end; % No scale
    
    if strcmp(order,'MP')
        mag=data_in(:,:,1);
        ph=data_in(:,:,2);
    elseif strcmp(order,'PM')
        mag=data_in(:,:,2);
        ph=data_in(:,:,1);
    end

    % Phase scaling parameter
    if isinteger(data_in)
        normlz=double(intmax(class(data_in)))-double(intmin(class(data_in)))+1;
        ph=single(ph); mag=single(mag);
    else
        normlz=1;
    end
    % Magnitude mapping
    switch mag_map
        case 'NS' % Nothing to do
        case 'LINM'
            mag = mag / lin_scale_factor;
        case 'LINP'
            mag = sqrt(mag / lin_scale_factor);
        case {'LOGM', 'LOGP'} % These seem identical from CMETAA spec.  10*log10(LINP)=10*log10(LINM^2)=20*log10(LINM)
            % The "DBperSTEP" value in the NITF CMETAA spec seems to be
            % required for magnitude mapping of the LOGM and LOGP mappings,
            % but it is not given in the CMETAA TRE, so it must be known
            % from elsewhere.  Here for LOGM and LOGP, we just use a
            % DBperSTEP that assumes data magnitude has been remapped to 8
            % bits from 16.
            db_per_step = 20 * log10(single(intmax('int16'))) / single(intmax('uint8'));
            mag = 10.^(mag * db_per_step / 20);
        case 'LLM'
            % For LLM, DBperSTEP must be set so that LINM equals LOGM at
            % the transition point.
            db_per_step = 20 * log10(ll_tp) / ll_tp;
            above_tp = mag>ll_tp;
            mag(above_tp) = 10.^(mag(above_tp) * db_per_step / 20); % LOGM portion
            mag = mag / lin_scale_factor; % LINM portion
        otherwise
            error('OPEN_CNITF_READER:UnsupportedFormat','Unsupported magnitude mapping.');
    end
    % Compute IQ
    data_out=mag.*exp(2*pi*1i*ph/normlz);
end

% LM C32 compression
function data_out = c32_to_complex(data_in)
% Lockheed-Martin C32 'compression' simply takes each 32-bit IEEE
% big-endian float in a complex number and lops off the last 16 bits of the
% mantissa. So since a typical 32-bit big-endian float looks like:
%
%     00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 ... 31
%     s  e  e  e  e  e  e  e  e  m  m  m  m  m  m  m  m  ... m
%
% where 's' indicates a sign bit, 'e' indicates an exponent bit, and 'm'
% indicates a mantissa bit (fraction bit) then a Lockheed-Martin C32 real
% or imaginary part is simply:
%
%     00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15    | Rest truncated !
%     s  e  e  e  e  e  e  e  e  m  m  m  m  m  m  m     |                !
%
% i.e. we are left with only 7 bits of real precision.  To rebuild this
% into a usable, if unreliable form, we simply insert unset (zero) bits to
% replace the missing ones:
%
%     00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 ... 31
%     s  e  e  e  e  e  e  e  e  m  m  m  m  m  m  m  0  ... 0

    data_out=zeros(size(data_in));
    data_out(:)=typecast(bitshift(uint32(data_in(:)),16),'single');
    data_out=complex(data_out(:,:,1),data_out(:,:,2));
end

% Determine the string describing the MATLAB class corresponding to the
% data stored in this image
function matlab_type_string = nitf2matlab_datatype(imagesubhdr)
    bitsperpixel=imagesubhdr.NBPP;
    switch imagesubhdr.PVTYPE
        case 'INT'
            typestring='uint';
        case 'SI'
            typestring='int';
        case 'R'
            typestring='float';
        case 'C'
            typestring='float';
            bitsperpixel=bitsperpixel/2;
        otherwise
            error('OPEN_CNITF_READER:UnsupportedFormat','Unsupported NITF format');
    end
            
    matlab_type_string=[typestring num2str(bitsperpixel)];
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////