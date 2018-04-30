function [ SICD_meta, NITF_meta ] = read_sicd_meta( filename )
%READ_SICD_META Read Sensor Independent Complex Data (SICD) metadata from 
%file into a MATLAB structure
%
% SICD (version >= 0.3) is stored in a NITF container.  NITF is a
% complicated format that involves lots of fields and configurations
% possibilities. Fortunately, SICD only really uses a small, specific
% portion of the NITF format.  This function extracts only the few parts of
% the NITF metadata necessary for reading a SICD NITF file.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Open file
NITF_meta = read_sicd_nitf_offsets(filename);
% SICD Volume 2, File Format Description, section 3.1.1 says that SICD XML
% metadata must be stored in first DES.
SICD_DES_Offset = NITF_meta.minimal.desOffsets(1);
SICD_DES_Length = NITF_meta.minimal.desLengths(1);

%% Read SICD XML metadata from the data extension segment
fid = fopen(filename);
fseek(fid,SICD_DES_Offset,'bof');
SICD_xml_string = fread(fid,SICD_DES_Length,'uint8=>char')';
fclose(fid);
SICD_meta = sicdxml2struct( xmlread( java.io.StringBufferInputStream( SICD_xml_string ) ) );

%% Adjust frequencies in metadata to be true, not offset values, if
% reference frequency is available
if isfield(SICD_meta,'RadarCollection') && isfield(SICD_meta.RadarCollection,'RefFreqIndex') && ...
        SICD_meta.RadarCollection.RefFreqIndex && exist('sicd_ref_freq','file')
    % Get reference frequency from function the user can place in MATLAB path
    SICD_meta = sicd_apply_ref_freq(SICD_meta, sicd_ref_freq());
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////