function  MENSRA  = readMENSRA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READMENSRA Create structure for metadata fields.
%   MENSRA = READMENSRA(FID) returns a structure of fields and their
%   associated values.

% Written by Matt Donath, NGA, matthew.b.donath@nga.ic.gov
% References:
%  > MIL-STD-2500A, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.0
%  > MIL-STD-2500C, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.1
%  > STDI-0001,     NATIONAL SUPPORT DATA EXTENSIONS (SDE)(VERSION 1.3/CN2)
%                   FOR THE NATIONAL IMAGERY TRANSMISSION FORMAT (NITF)
%  > STDI-0002,     THE COMPENDIUM OF CONTROLLED EXTENSIONS (CE) FOR THE 
%                   NATIONAL IMAGERY TRANSMISSION FORMAT (NITF) VERSION 2.1
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Create the metadata structure for MENSRA.

MENSRA.CETAG = fread(fid,6,'uint8=>char')';
MENSRA.CEL = str2double(fread(fid,5,'uint8=>char')');

MENSRA.CP_LOC = fread(fid,21,'uint8=>char')';
MENSRA.CP_ALT = str2double(fread(fid,6,'uint8=>char')');
MENSRA.OF_PC_R = str2double(fread(fid,7,'uint8=>char')');
MENSRA.OF_PC_A = str2double(fread(fid,7,'uint8=>char')');
MENSRA.COSGRZ = str2double(fread(fid,7,'uint8=>char')');
MENSRA.RGCCRP = str2double(fread(fid,7,'uint8=>char')');
MENSRA.RLMAP = fread(fid,1,'uint8=>char')';
MENSRA.CCRP_ROW = str2double(fread(fid,5,'uint8=>char')');
MENSRA.CCRP_COL = str2double(fread(fid,5,'uint8=>char')');
MENSRA.ACFT_LOC = fread(fid,21,'uint8=>char')';
MENSRA.ACFT_ALT = str2double(fread(fid,5,'uint8=>char')');
MENSRA.C_R_NC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_R_EC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_R_DC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_AZ_NC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_AZ_EC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_AZ_DC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_AL_NC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_AL_EC = str2double(fread(fid,7,'uint8=>char')');
MENSRA.C_AL_DC = str2double(fread(fid,7,'uint8=>char')');

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////