function [ readerobj ] = open_cos_reader( filename )
%OPEN_COS_READER Intiates a reader object for TerraSAR-X COSAR file format.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


%% Defaults if no XML found
symmetry=[0 0 1]; % COSAR written in range lines
meta=struct();

%% Find XML metadata if available
[fullpath, fileonly, ext]=fileparts(filename);
[parentpath,pathname]=fileparts(fullpath);
if strcmp(pathname,'IMAGEDATA')
    xmlfile=dir(fullfile(parentpath,'*.xml'));
    geoxmlfile=dir(fullfile(parentpath,'ANNOTATION','GEOREF.xml'));
    if (length(xmlfile)==1)&&istsx(fullfile(parentpath,xmlfile.name))
        tsxxml_meta=xmlread(fullfile(parentpath,xmlfile.name));
        if ~isempty(geoxmlfile)
            tsxgeoxml_meta = xmlread(fullfile(parentpath,'ANNOTATION','GEOREF.xml'));
        else
            tsxgeoxml_meta = [];
        end
        meta=meta2sicd_tsxxml(tsxxml_meta, tsxgeoxml_meta);
        % XML describes all images of a polarimetric collection
        % Select metadata associated with only this COSAR file
        xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
        layerIndex=...
            str2double(xp.evaluate(...
            ['level1Product/productComponents/imageData[file/location/filename="' fileonly ext '"]/@layerIndex'],...
            tsxxml_meta));
        meta=meta{layerIndex};
        meta.native.tsxxml=tsxxml_meta;
        symmetry(1) = (meta.SCPCOA.SideOfTrack=='L'); % Flip for left looking
        % Reading ValidData region from COSAR is slow, so we comment it out
        % meta.ImageData.ValidData = cosar_valid_data(filename,symmetry(1));
        % valid_latlons=point_slant_to_ground(...
        %     [[meta.ImageData.ValidData.Vertex.Row]; ...
        %     [meta.ImageData.ValidData.Vertex.Col]], meta);
        % for j = 1:size(valid_latlons,2)
        %     meta.GeoData.ValidData.Vertex(j).Lat=valid_latlons(1,j);
        %     meta.GeoData.ValidData.Vertex(j).Lon=valid_latlons(2,j);
        % end
    end
end

%% Open reader
readerobj=open_cos_reader_noxml(filename,symmetry);
meta = setstructfields(meta,readerobj.get_meta()); % Merge COSAR and XML metadata
readerobj.get_meta=@() meta;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////