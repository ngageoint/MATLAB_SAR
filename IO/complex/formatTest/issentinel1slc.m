function boolout = issentinel1slc( filename )
%ISSENTINEL1SLC Test for Sentinel-1 SLC "SAFE" file
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Might be cell array of SAFE file ane orbit file
if iscell(filename)
    boolout = issentinel1slc(filename{1}); % Check validity of SAFE file
    % Check for validity of orbit file here
    if boolout && mightbexml(filename{2})
        xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
        domnode = xmlread(filename{2});
        try
            % The local-name XPath construct is used to avoid dealing with
            % fields that use namespaces.
            boolout = ...
                strncmpi(char(xp.evaluate(... % Is a Sentinel-1 EOF file
                    'Earth_Explorer_File/Earth_Explorer_Header/Fixed_Header/Mission',domnode)),'SENTINEL-1', 10) && ...
                any(strcmp(char(xp.evaluate(... % Is a compatible orbit file type
                    'Earth_Explorer_File/Earth_Explorer_Header/Fixed_Header/File_Type',domnode)),{'AUX_RESORB','AUX_POEORB'}));
        catch
            boolout=false;
        end
    end
    return;
end

% Or just a SAFE file
if mightbexml(filename)
    xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
    domnode = xmlread(filename);
    % Overall wrapping node is in a namespace that makes things messier
    % (than they already are) in a XPath call, so we go down one level.
    domnode = domnode.getFirstChild;
    try
        % The local-name XPath construct is used to avoid dealing with
        % fields that use namespaces.
        boolout = ...
            strcmp(char(xp.evaluate(... % Is a Sentinel-1 SAFE file
                ['metadataSection/' ...
                'metadataObject[@ID="platform"]/' ...
                'metadataWrap/' ...
                'xmlData/' ...
                '*[local-name()=''platform'']/' ...
                '*[local-name()=''familyName'']'],domnode)),'SENTINEL-1') && ...
            strcmp(char(xp.evaluate(... % Is an SLC product type
                ['metadataSection/' ...
                'metadataObject[@ID="generalProductInformation"]/' ...
                'metadataWrap/' ...
                'xmlData/' ...
                '*[local-name()=''standAloneProductInformation'']/' ...
                '*[local-name()=''productType'']'],domnode)),'SLC');
    catch
        boolout=false;
    end
else
    boolout=false;
end

end


% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////