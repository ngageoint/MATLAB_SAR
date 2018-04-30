function [ readerobj ] = open_rs_reader( filename )
%OPEN_RS_READER Intiates a reader object for RADARSAT-2 or RCM product.xml file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Read and process XML metadata
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
domnode=xmlread(filename);
symmetry=[0 0 0];
% Because RS2/RCM XML uses a default namespace, we have to use a local-name
% query, rather than the much simpler traditional XPath query that assumes
% no namespace.  See isrs2.m for more details.
collector_name=char(xp.evaluate(['/*[local-name()=''product'']'...
    '/*[local-name()=''sourceAttributes'']'...
    '/*[local-name()=''satellite'']'], domnode));
if strcmpi(collector_name, 'RADARSAT-2')
    ia_str = 'imageAttributes';
    ipdf_str=['/*[local-name()=''product'']'...
        '/*[local-name()=''imageAttributes'']'...
        '/*[local-name()=''fullResolutionImageData'']'];
    beta_lut_str=['/*[local-name()=''product'']'...
        '/*[local-name()=''imageAttributes'']'...
        '/*[local-name()=''lookupTable''][@incidenceAngleCorrection="Beta Nought"]'];
    beta_lut_str = fullfile(fileparts(filename), char(xp.evaluate(beta_lut_str, domnode)));
    noise_str = ''; % No separate noise file for RS2
elseif strncmpi(collector_name, 'RCM', 3)
    ia_str = 'imageReferenceAttributes';
    ipdf_str=['/*[local-name()=''product'']'...
        '/*[local-name()=''sceneAttributes'']'...
        '/*[local-name()=''imageAttributes'']'...
        '/*[local-name()=''ipdf'']'];
    beta_lut_str=['/*[local-name()=''product'']'...
        '/*[local-name()=''imageReferenceAttributes'']'...
        '/*[local-name()=''lookupTableFileName''][@sarCalibrationType="Beta Nought"]'];
    beta_lut_str = fullfile(fileparts(filename), 'calibration', char(xp.evaluate(beta_lut_str, domnode)));
    noise_str=['/*[local-name()=''product'']'...
        '/*[local-name()=''imageReferenceAttributes'']'...
        '/*[local-name()=''noiseLevelFileName'']'];
    noise_str = fullfile(fileparts(filename), 'calibration', char(xp.evaluate(noise_str, domnode)));
end
lineOrder=char(xp.evaluate(['/*[local-name()=''product'']'...
    '/*[local-name()=''' ia_str ''']'...
    '/*[local-name()=''rasterAttributes'']'...
    '/*[local-name()=''lineTimeOrdering'']'],...
    domnode));
lookdir=upper(char(xp.evaluate(['/*[local-name()=''product'']'...
    '/*[local-name()=''sourceAttributes'']'...
    '/*[local-name()=''radarParameters'']'...
    '/*[local-name()=''antennaPointing'']'],...
    domnode)));
sampleOrder=char(xp.evaluate(['/*[local-name()=''product'']'...
    '/*[local-name()=''' ia_str ''']'...
    '/*[local-name()=''rasterAttributes'']'...
    '/*[local-name()=''pixelTimeOrdering'']'],...
    domnode));
symmetry(1)=xor(strcmp(lineOrder,'Decreasing'),lookdir(1)=='L');
symmetry(2)=strcmp(sampleOrder,'Decreasing');
if exist(beta_lut_str, 'file')
    betanode={xmlread(beta_lut_str)};
else
    betanode = {};
end
if exist(noise_str, 'file')
    try % Simulated RCM datasets have misformed noise XML
        betanode = {betanode{:} xmlread(noise_str)};
    end
end
meta=meta2sicd_rs_xml(domnode,betanode{:});

%% Get all files from a polarimetric collect
num_files=str2double(xp.evaluate(['count(' ipdf_str ')'],domnode));
readerobj=cell(1,num_files);
for i=1:num_files
    basepathname=fileparts(filename);
    datafilename=char(xp.evaluate([ipdf_str '[' num2str(i) ']'],domnode));
    % If format is needed, do this:
    % format=char(xp.evaluate(['/*[local-name()=''product'']'...
    %                          '/*[local-name()=''' ia_str ''']'...
    %                          '/*[local-name()=''productFormat'']']),...
    %                         domnode);
    % Return individual reader objects for each polarimetric channel
    readerobj{i}=open_tiff_reader_noxml(fullfile(basepathname,datafilename),symmetry);
    meta{i}.native.tiff=readerobj{i}.get_meta();
    readerobj{i}.get_meta = @() meta{i};
end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////