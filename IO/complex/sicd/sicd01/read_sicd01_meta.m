function [ meta ] = read_sicd01_meta( filename )
%READ_SICD01_META Read metadata from Sensor Independent Complex Data (SICD)
%file, version 0.1
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open file
fid=fopen(filename,'r','b');
if(fid<0)
    error('READ_SICD01_META:InvalidFile','Unable to open file.');
end
% Read in non-xml header data
meta.format=fgetl(fid);
if(~strcmp(meta.format,'SICD-E/0.1'))
    error('READ_SICD01_META:UnrecognizedVersion','Unrecognized version of SICD format.  Only version 0.1 supported.');
end
for i=1:10
    [fieldname,value]=strtok(fgetl(fid),':=');
    fieldname=regexprep(fieldname,'-','_'); % Make sure field is valid Matlab variable name
    meta.header.(fieldname)=value(3:end);
    doubval=str2double(meta.header.(fieldname));
    if(~isnan(doubval)) % If value is numeric, store as such
        meta.header.(fieldname)=doubval;
    end
end
meta.header_length_in_bytes=ftell(fid);
fread(fid,2,'uchar'); % Blank line (/12/10) in between header and xml
% Read XML
SICD_xml_string = fread(fid,meta.header.xml_content_length,'uchar=>char',0)';
meta.xmldata = xml2simplestruct( xmlread( java.io.StringBufferInputStream( SICD_xml_string ) ) );
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////