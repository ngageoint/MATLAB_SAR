function boolout = istsx( filename )
%ISTSX Test for TerrSAR-X XML description file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if mightbexml(filename)
    xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();
    try
        boolout=isequal(regexp(char(xp.evaluate('level1Product/generalHeader/mission',xmlread(filename))),'T[SD]X-1'),1);
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