function write_sicd_dessubhdr(obj)
%WRITE_SICD_DESSUBHDR Writes the NITF Data Extension Segment sub-header for a SICD
%
% It's not expected that the typical user would call this function, rather
% it is called as part of a bigger SICD writer.
%
% Inputs:
%       obj:       This is the 'hidden' object-specific handle.  It is
%                  similar to 'this' in Java.
%
% References:
%       MIL-STD-2500C, Department of Defense Interface Standard
%                      National Imagery Transmission Format version 2.1
%                      01 May 2006
%
% Written by: Tom Krauss and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Always writing version 2.1 NITF

fwriten(obj.FID, 'DE', 2);   % DE
fwriten(obj.FID, 'XML_DATA_CONTENT', 25);  % DESID
fwriten(obj.FID, '01', 2);   % DESVER
obj.write_sicd_security_tags(); % Security tags go here
%fwriten(obj.FID, '0000', 4);   % DESSHL - Complete inclusion of all User-defined Subheader Subfield was not required prior to SICD version 1.0
fwriten(obj.FID, '0773', 4);   % DESSHL
fwriten(obj.FID, '99999', 5); % DESCRC - CRC not computed
fwriten(obj.FID, 'XML', 8); % DESSHFT
if isfield(obj.sicdmeta,'ImageCreation') && isfield(obj.sicdmeta.ImageCreation, 'DateTime')
    fwriten(obj.FID, datestr(obj.sicdmeta.ImageCreation.DateTime,'yyyy-mm-ddTHH:MM:SSZ'), 20); % DESSHDT
else
    fwriten(obj.FID, '', 20); % DESSHDT
end
fwriten(obj.FID, '', 40); % DESSHRP
fwriten(obj.FID, 'SICD Volume 1 Design & Implementation Description Document', 60); % DESSHSI
if isfield(obj.sicdmeta,'schema') && isfield(obj.sicdmeta.schema,'version')
    schema_version = obj.sicdmeta.schema.version;
else
    schema_version = '1.3.0';  % Default to current
end
if isfield(obj.sicdmeta,'schema') && isfield(obj.sicdmeta.schema,'date')
    schema_date = obj.sicdmeta.schema.date;
else
    schema_date = '2021-11-30T00:00:00Z';  % Default to current
end
if isfield(obj.sicdmeta,'schema') && isfield(obj.sicdmeta.schema,'urn')
    schema_urn = obj.sicdmeta.schema.urn;
else
    schema_urn = 'urn:SICD:1.3.0';  % Default to current
end
fwriten(obj.FID, schema_version, 10); % DESSHSV
fwriten(obj.FID, schema_date, 20); % DESSHSD
fwriten(obj.FID, schema_urn, 120); % DESSHTN
if isfield(obj.sicdmeta,'GeoData') && isfield(obj.sicdmeta.GeoData,'ImageCorners') && ...
    isfield(obj.sicdmeta.GeoData.ImageCorners,'ICP')
    DESSHLPG = sprintf('%+012.8f%+013.8f%+012.8f%+013.8f%+012.8f%+013.8f%+012.8f%+013.8f',...
        obj.sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lat,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lon,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lat,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.FRLC.Lon,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lat,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.LRLC.Lon,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lat,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.LRFC.Lon,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lat,...
        obj.sicdmeta.GeoData.ImageCorners.ICP.FRFC.Lon);
else
    DESSHLPG = '';
end
fwriten(obj.FID, DESSHLPG, 125); % DESSHLPG
fwriten(obj.FID, '', 25); % DESSHLPT
fwriten(obj.FID, '', 20); % DESSHLI
fwriten(obj.FID, '', 120); % DESSHLIN
fwriten(obj.FID, '', 200); % DESSHABS

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////