function boolout = isrcmmanifest( filename )
%ISRCM Test for Radar Constellation Mission "SAFE" manifest file
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Check for valid SAFE file
if mightbexml(filename)
    xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
    domnode = xmlread(filename);
    % Overall wrapping node is in a namespace that makes things messier
    % (than they already are) in a XPath call, so we go down one level.
    domnode = domnode.getFirstChild;
    try
        % The local-name XPath construct is used to avoid dealing with
        % fields that use namespaces.
        productfilename = fullfile(fileparts(filename),'metadata','product.xml');
        boolout = (strcmp(char(xp.evaluate(... % Is a Sentinel-1 SAFE file
                ['informationPackageMap/' ...
                '*[local-name()=''contentUnit'']/' ...
                '@unitType'],domnode)),'RCM Product Information Package') && ...
            strcmp(char(xp.evaluate(... % Is an SLC product type
                ['informationPackageMap/' ...
                '*[local-name()=''contentUnit'']/' ...
                '@textInfo'],domnode)),'RCM SLC Product') && ...
            exist(productfilename, 'file') && ...
            isrs(productfilename));
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