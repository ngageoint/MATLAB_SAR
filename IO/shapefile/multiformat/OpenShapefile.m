function fid = OpenShapefile(filename,format,varargin)
%OPENSHAPEFILE Writes header to specified shapefile format 
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Function writes header to Socet GRFX .grfx file, RemoteView .rvl file or
% Google Earth .kml file
%
%CLASSIFICATION: 
%   Unclassified
%
%INPUTS:
%  filename      - required: shapefile name
%  format        - required: 'grfx','rvl' or 'kml'
%  PropertyName  - The PropertyName and PropertyValue are key-value
%  PropertyValue - pairs used throughout the tool kit...
%
% Allowed properties:
%       Property name         Description
%       -----------------------------------------------------------------
%       name            The name to put into the KML file.  This name is
%                       displayed in the "Places" tree viewer by Google
%                       Earth.
%
%OUTPUTS:
%  fid - open file handle to shapefile
% 
%VERSION:
%   1.0 
%     - Tim Cox 20090324
%     - initial version


p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('name', 'Untitled');
p.FunctionName = 'OPENSHAPEFILE';
p.parse(varargin{:});

fid = fopen(filename,'w');

format = lower(format);

if (strcmp(format,'grfx') == 1)
    fprintf(fid,'<?xml version="1.0" ?>\n');
    fprintf(fid,'<Graphics appVersionMajor="2" appVersionMinor="3" appVersionBuild="1"><GmGraphicsLayout m_has_changed="true" m_manage_memory="true">\n');
    fprintf(fid,'<GmLayerList name="m_layer_list">\n');
    fprintf(fid,'<GmLayer m_layer_name="GRFXLayer" m_manage_memory="true" m_refresh_indices="true" m_exclusive="false">\n');
    fprintf(fid,'<GrGmGraphicList name="m_graphic_list">\n\n');
    return;
elseif (strcmp(format,'rvl') == 1)
    fprintf(fid,'<<\n');
    fprintf(fid,'Layer	<<\n');
    fprintf(fid,'DisplayList	[\n\n');
    return;
elseif (strcmp(format,'kml') == 1)
    fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
	  fprintf(fid,'<kml xmlns = "http://www.opengis.net/kml/2.2">\n');
	  fprintf(fid,'<Document>\n\n');
	  fprintf(fid,'<name>%s</name>\n\n', p.Results.name);
    return;
else    
    return;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

