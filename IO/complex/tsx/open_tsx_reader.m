function [ readerobj ] = open_tsx_reader( filename )
%OPEN_TSX_READER Intiates a reader object for TSX XML description file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Read in XML metadata
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
domnode=xmlread(filename); % Read main TSX XML descriptor file
filepath=fileparts(filename);
geoxmlfile=dir(fullfile(filepath,'ANNOTATION','GEOREF.xml'));
if ~isempty(geoxmlfile)
    geo_dom_node = xmlread(fullfile(filepath,'ANNOTATION','GEOREF.xml'));
else
    geo_dom_node = [];
end
meta=meta2sicd_tsxxml(domnode, geo_dom_node);
for i=1:length(meta)
    meta{i}.native.tsxxml=domnode;
end
symmetry=[(meta{1}.SCPCOA.SideOfTrack=='L') 0 1]; % COSAR written in range lines

%% Setup readers for each component
num_files=str2double(xp.evaluate('count(level1Product/productComponents/imageData/file/location)',domnode));
filelist=cell(1,num_files); readerobj=filelist; % Create empty cell arrays
for i=1:num_files
    basepathname=fileparts(filename);
    xpath_str_base=['level1Product/productComponents/imageData[' num2str(i) ']/file/'];
    host=char(xp.evaluate([xpath_str_base 'location/host'],domnode));
    pathname=char(xp.evaluate([xpath_str_base 'location/path'],domnode));
    datafilename=char(xp.evaluate([xpath_str_base 'location/filename'],domnode));
    filelist{i}=fullfile(basepathname,host,pathname,datafilename);
    format = char(xp.evaluate('level1Product/productInfo/imageDataInfo/imageDataFormat',domnode));
    if strcmp(format,'COSAR')
        readerobj{i}=open_cos_reader_noxml(filelist{i},symmetry);
        cos_meta=readerobj{i}.get_meta();
        meta{i}.native.cos=cos_meta.native.cos;
        if (meta{i}.ImageData.NumCols~=cos_meta.native.cos.az)||...
            (meta{i}.ImageData.NumRows~=cos_meta.native.cos.rs)
            error('OPEN_TSX_READER:INCONSISTENT_SIZE','COSAR and XML metadata do not match.');
        end
        % Reading ValidData region from COSAR is slow, so we comment it out
        % meta{i}.ImageData.ValidData = cosar_valid_data(filelist{i},symmetry(1));
        % valid_latlons=point_slant_to_ground(...
        %     [[meta{i}.ImageData.ValidData.Vertex.Row]; ...
        %     [meta{i}.ImageData.ValidData.Vertex.Col]], meta{i});
        % for j = 1:size(valid_latlons,2)
        %     meta{i}.GeoData.ValidData.Vertex(j).Lat=valid_latlons(1,j);
        %     meta{i}.GeoData.ValidData.Vertex(j).Lon=valid_latlons(2,j);
        % end
    elseif strcmp(format,'GEOTIFF')
        readerobj{i}=open_tiff_reader(filelist{i}); % Not necessarily oriented properly
        meta{i}.native.tiff=readerobj{i}.get_meta();
    end
    readerobj{i}.get_meta = @() meta{i};
end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////