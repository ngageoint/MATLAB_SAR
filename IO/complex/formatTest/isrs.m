function boolout = isrs( filename )
%ISRS Test for RADARSAT-2 or RCM product.xml description file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if mightbexml(filename)
    xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
    try
        satellite=xp.evaluate(['/*[local-name()=''product'']'...
            '/*[local-name()=''sourceAttributes'']'...
            '/*[local-name()=''satellite'']'],xmlread(filename));
        % OK, now a long explanation for a single line of code.  The above
        % XPath expression is less efficient (and much more clumsy to
        % write) than the following which worked in MATLAB 2009b and below:
        % boolout=strcmpi('RADARSAT-2',...
        %     xp.evaluate('product/sourceAttributes/satellite',xmlread(filename)));
        % However, XPath with no namespace resolver merely queries for an
        % element that is not associated with any namespace.  Since RS2 XML
        % files have a default namespace, every element is associated with
        % a namespace, and a basic XPath query like that will not work.
        % Actually the basic XPath query did work in MATLAB 2009b and
        % before because the Sun XPath parser used assumed the default
        % namespace, if none was given.  In MATLAB 2010a, Mathworks added
        % the Saxon XPath parser which actually behaves closer to XPath
        % spec, in my understanding, and the simple solution no longer
        % worked.
        boolout=strcmpi('RADARSAT-2', satellite) || ...
            strncmpi('RCM', satellite, 3);
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