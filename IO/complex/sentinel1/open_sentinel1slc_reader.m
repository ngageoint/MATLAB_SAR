function [ readerobj ] = open_sentinel1slc_reader( filename )
%OPEN_SENTINEL1SLC_READER Intiates a reader object for Sentinel-1 SAFE file
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Constants
symmetry=[0 0 0]; % Sentinel-1 is only right-looking, so this should be true for all datasets
% The above symmetry is true for MATLAB's TIFF reader.  However, MATLAB's
% TIFF reader won't read the Sentinel-1 complex TIFFs, so we had to make
% our own, which requires an adjustment to the symmetry, since our homemade
% TIFF reader returns data transposed from MATLAB's (because that's the
% order of the data in the file and how it should be read.)
symmetry(3) = ~symmetry(3);
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

%% Allow for SAFE and EOF files to be passed in cell array
if iscell(filename)
    eof_filename = filename{2};
    filename = filename{1};
end

%% Read and process XML metadata
manifest_domnode = xmlread(filename);
files = s1manifest_files(manifest_domnode);
meta_manifest = meta2sicd_s1manifest(manifest_domnode);

%% Create reader objects for all measurement files
% We should create a readerobject for each burst, of which there may be
% several in each measurement file.
readerobj = {};
meta = {};
for i=1:numel(files) % There will be a file for each swath and polarization
    basepathname = fileparts(filename);
    datafilename = files(i).data;
    meta_tiff = read_tiff_tags(fullfile(basepathname,datafilename));
    if ~isempty(files(i).product) && ...
            exist(fullfile(basepathname, files(i).product),'file')
        product_domnode = xmlread(fullfile(basepathname, files(i).product));
        % meta2sicd_s1product will return an array of metadata structures if there are multiple bursts
        if exist('eof_filename','var')
            eof_domnode = xmlread(eof_filename);
            meta_product = meta2sicd_s1product(product_domnode, eof_domnode);
        else
            meta_product = meta2sicd_s1product(product_domnode);
        end
        if ~isempty(files(i).calibration) && ...
                exist(fullfile(basepathname, files(i).calibration),'file')
            meta_cal = meta2sicd_s1cal(xmlread(fullfile(basepathname, files(i).calibration)));
            if isfield(meta_cal,'Radiometric')
                for j = 1:numel(meta_product)
                    meta_product{j}.Radiometric = meta_cal.Radiometric;
                    % meta2sicd_s1cal only populates BetaZeroSFPoly.  Next
                    % we derive the rest.
                    meta_product{j} = derived_sicd_fields(meta_product{j});
                end
            end
        end
        if ~isempty(files(i).noise) && ...
                exist(fullfile(basepathname, files(i).noise),'file')
            meta_noise = meta2sicd_s1noise(xmlread(fullfile(...
                basepathname, files(i).noise)), meta_product);
            for j = 1:numel(meta_product)
                if iscell(meta_noise) && ...
                        isfield(meta_noise{j},'Radiometric') && ...
                        isfield(meta_noise{j}.Radiometric,'NoiseLevel')
                    meta_product{j}.Radiometric.NoiseLevel = ...
                        meta_noise{j}.Radiometric.NoiseLevel;
                end
            end
        end
        num_bursts = numel(meta_product);
        if num_bursts == 1 % Stripmap, single burst, open entire file
            readerobj{i} = open_ctiff_reader_noxml(fullfile(basepathname,datafilename),symmetry);
        else % Interferometric wide, multiple swaths (files), and multiple bursts within a file
            base_reader = open_ctiff_reader_noxml(...
                    fullfile(basepathname,datafilename),symmetry);
            num_lines_burst = str2double(char(xp.evaluate(...
                'product/swathTiming/linesPerBurst',product_domnode)));
            for j = 1:num_bursts
                readerobj{end+1} = subset_reader(base_reader,...
                    num_lines_burst*(j-1) + [1 num_lines_burst],[1 meta_tiff{1}.ImageWidth]);
            end
        end
        for j = 1:num_bursts
            % Handle dual-polarization case.  Label channel number
            % appropriately for ordering in manifest file.
            meta_product{j}.ImageFormation.RcvChanProc.ChanIndex = find(strcmp(...
                {meta_manifest.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization},...
                meta_product{j}.ImageFormation.TxRcvPolarizationProc));
            meta{end+1} = setstructfields(meta_manifest, meta_product{j});
            % meta should already be set to this from meta_product:
            % meta{end}.ImageData.NumRows = meta_tiff{1}.ImageWidth;
            meta{end}.native.tiff = meta_tiff;
            readerobj{numel(meta)}.get_meta = @() meta{end};
        end
    else % No annotation metadata could be found
        meta{end+1} = meta_manifest;
        meta{end}.ImageData.NumCols = meta_tiff{1}.ImageLength;
        meta{end}.ImageData.NumRows = meta_tiff{1}.ImageWidth;
        meta{end}.native.tiff = meta_tiff;
        readerobj{i}.get_meta = @() meta{end};
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////