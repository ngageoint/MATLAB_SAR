function boolout = issicd(filename)
% Sensor independent complex data (version >= 0.3)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

try
    nitf_meta = read_sicd_nitf_offsets(filename);
    fid = fopen(filename,'r');
    fseek(fid, nitf_meta.minimal.desOffsets(1), 'bof'); % First DES must be SICD XML
    domnode = xmlread( java.io.StringBufferInputStream( fread(fid,nitf_meta.minimal.desLengths(1),'uint8=>char')));
    fclose(fid);
    boolout = strncmpi('SICD', domnode.getDocumentElement.getNodeName,4);
catch
    % read_sicd_nitf_offsets throws an error for non-NITF files
    % non-XML DES will also throw error
    boolout = false; % Not SICD either way
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////