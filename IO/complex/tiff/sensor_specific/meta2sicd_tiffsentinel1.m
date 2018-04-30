function [ meta, symmetry ] = meta2sicd_tiffsentinel1( tiff_filename )
%META2SICD_TIFFSENTINEL1 Compile SICD metadata for Sentinel-1 data package
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

symmetry=[0 0 0]; % Data is time-ordered, but always right-looking, so this is equivalent to view-from-above

pathname = fileparts(tiff_filename);
s1xmlfile=dir(fullfile(fileparts(pathname),'manifest.safe')); % Sentinel 1 file
if (length(s1xmlfile)==1)&&issentinel1slc(fullfile(fileparts(pathname),s1xmlfile.name))
    % Process manifest.safe file
    manifest_path = fileparts(pathname);
    manifest_domnode=xmlread(fullfile(manifest_path,'manifest.safe'));
    meta_manifest=meta2sicd_s1manifest(manifest_domnode);
    % Which measurement file in the manifest did we open?
    files = s1manifest_files(manifest_domnode);
    index = strcmp(fullfile(manifest_path,{files.data}),tiff_filename);
    if ~isempty(files(index).product) && ...
            exist(fullfile(manifest_path, files(index).product),'file')
        % Find associated xml metadata file
        product_domnode=xmlread(fullfile(manifest_path, files(index).product));
        meta_product = meta2sicd_s1product(product_domnode);
        for j = 1:numel(meta_product)
            % Handle dual-polarization case.  Label channel number appropriate
            % for ordering in manifest file.
            meta_product{j}.ImageFormation.RcvChanProc.ChanIndex = find(strcmp(...
                {meta_manifest.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization},...
                meta_product{j}.ImageFormation.TxRcvPolarizationProc));
            meta{j} = setstructfields(meta_manifest, meta_product{j});
        end
    else % Unable to find annotation file
        meta = meta_manifest;
    end
else
    meta=struct();
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////