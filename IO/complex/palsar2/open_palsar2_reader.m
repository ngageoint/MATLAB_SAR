function [ readerobj ] = open_palsar2_reader( filename )
%OPEN_PALSAR2_READER Intiates a reader object for ALOS PALSAR 2 dataset
%
% Handles ALOS PALSAR 2 level 1.1 (SLC) data
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup reader
filenames = palsar2_files(filename);
if any(strcmp(filename, filenames.IMG))
    % If a single polarization is passed, only open that polarization.  If
    % a VOL/LED/TRL file is passed, open all polarizations.
    filenames.IMG = {filename};
end
meta_base = struct();
if isfield(filenames,'VOL')
    meta_base.native.vol = read_ceos_vol_meta(filenames.VOL);
    meta_base = setstructfields(meta_base, meta2sicd_palsar2vol(meta_base.native.vol));
end
if isfield(filenames,'LED')
    meta_base.native.led = read_ceos_led_meta(filenames.LED);
    meta_base = setstructfields(meta_base, meta2sicd_palsar2led(meta_base.native.led));
end
% Not really anything useful for us in the TRL file, so we skip that one
for i = 1:numel(filenames.IMG)
    img_meta = read_ceos_img_meta(filenames.IMG{i});
    datasize{i}=[img_meta.num_pixels img_meta.num_lines];
    data_offset{i}=[img_meta.rec_length + img_meta.prefix_bytes... % Byte offset to first data sample 
                 img_meta.suffix_bytes + img_meta.prefix_bytes]; % Byte spacing between rows
    meta{i} = setstructfields(meta_base, meta2sicd_palsar2img(img_meta));
    % Much of the useful metadata requires both he LED and IMG files to compute.
    if isfield(filenames,'LED')
        meta{i} = setstructfields(meta{i}, ...
            meta2sicd_palsar2ledimg(meta_base.native.led, img_meta));
        % Now that sensor model fields have been populated, we can populate
        % GeoData.SCP more precisely.
        ecf = point_image_to_ground([meta{i}.ImageData.SCPPixel.Row;meta{i}.ImageData.SCPPixel.Col],meta{i});
        meta{i}.GeoData.SCP.ECF.X=ecf(1);
        meta{i}.GeoData.SCP.ECF.Y=ecf(2);
        meta{i}.GeoData.SCP.ECF.Z=ecf(3);
        llh=ecf_to_geodetic([meta{i}.GeoData.SCP.ECF.X meta{i}.GeoData.SCP.ECF.Y meta{i}.GeoData.SCP.ECF.Z]);
        meta{i}.GeoData.SCP.LLH.Lat=llh(1);
        meta{i}.GeoData.SCP.LLH.Lon=llh(2);
        meta{i}.GeoData.SCP.LLH.HAE=llh(3);
    end
    meta{i} = derived_sicd_fields(meta{i});
    meta{i}.native.img = img_meta;
end
% Consolidate polarizations across images
pols = cellfun(@(x) x.ImageFormation.TxRcvPolarizationProc, meta, 'UniformOutput', false);
tpols = unique(cellfun(@(x) x(1), pols));
for i = 1:numel(filenames.IMG)
    meta{i}.ImageFormation.RcvChanProc.ChanIndex = i;
    if isscalar(tpols)
        meta{i}.RadarCollection.TxPolarization = tpols;
    else
        meta{i}.RadarCollection.TxPolarization = 'SEQUENCE';
        for j = 1:numel(tpols)
            meta{i}.RadarCollection.TxSequence.TxStep(j).TxPolarization = tpols(j);
        end
    end
    for j = 1:numel(pols)
        meta{i}.RadarCollection.RcvChannels(j).ChanParameters.TxRcvPolarization = pols{j};
    end
end
if isfield(meta{1}, 'SCPCOA') && isfield(meta{1}.SCPCOA, 'SideOfTrack') && ...
        upper(meta{1}.SCPCOA.SideOfTrack(1))=='L'
    symmetry=[0 1 1]; % PALSAR written in range lines
else
    symmetry=[0 0 1]; % PALSAR written in range lines
end
if strncmp(img_meta.sar_datatype_code, 'C*8', 3) % Assumes all IMG files are of same type
    datatype='float32';
    complextype=true;
else
    datatype='uint16';
    complextype=false;
end
endian='b';
bands=1;

%% Open reader
for i = 1:numel(filenames.IMG)
    readerobj{i}=open_generic_reader(filenames.IMG{i}, datasize{i},...
        datatype, complextype, data_offset{i}, endian, symmetry, bands);
    readerobj{i}.get_meta=@() meta{i};
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////