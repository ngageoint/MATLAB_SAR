function CloseShapefile(fid,format)
%CLOSESHAPEFILE Writes footer to specified shapefile format 
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Function writes footer to Socet GRFX .grfx file, RemoteView .rvl file or
% Google Earth .kml file
%
%CLASSIFICATION: 
%   Unclassified
%
%INPUTS:
%  fid      - required: File Handle of open Shape File
%  format   - required: 'grfx','rvl' or 'kml'
%
%OUTPUTS:
% 
%VERSION:
%   1.0 
%     - Tim Cox 20090324
%     - initial version

format = lower(format);

if (strcmp(format,'grfx') == 1)
    fprintf(fid,'</GrGmGraphicList>\n');
	fprintf(fid,'</GmLayer>\n');
	fprintf(fid,'</GmLayerList>\n');
	fprintf(fid,'</GmGraphicsLayout>\n');
	fprintf(fid,'</Graphics>\n');
    fclose(fid);
    return;
elseif (strcmp(format,'rvl') == 1)
    fprintf(fid,']\n');
	fprintf(fid,'>>\n');
	fprintf(fid,'>>\n');
    fclose(fid);
    return;
elseif (strcmp(format,'kml') == 1)
    fprintf(fid,'</Document>\n');
	fprintf(fid,'</kml>\n');
    fclose(fid);
    return;
else    
    return;
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
