function  MENSRB  = readMENSRB(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READMENSRB Create structure for metadata fields.
%   MENSRB = READMENSRB(FID) returns a structure of fields and their
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

%% Create the metadata structure for MENSRB.

MENSRB.CETAG = fread(fid,6,'uint8=>char')';
MENSRB.CEL = str2double(fread(fid,5,'uint8=>char')');

MENSRB.ACFT_LOC = fread(fid,25,'uint8=>char')';
MENSRB.ACFT_LOC_ACCY = str2double(fread(fid,6,'uint8=>char')');
MENSRB.ACFT_ALT = str2double(fread(fid,6,'uint8=>char')');
MENSRB.RP_LOC = fread(fid,25,'uint8=>char')';
MENSRB.RP_LOC_ACCY = str2double(fread(fid,6,'uint8=>char')');
MENSRB.RP_ELV = str2double(fread(fid,6,'uint8=>char')');
MENSRB.OF_PC_R = str2double(fread(fid,7,'uint8=>char')');
MENSRB.OF_PC_A = str2double(fread(fid,7,'uint8=>char')');
MENSRB.COSGRZ = str2double(fread(fid,7,'uint8=>char')');
MENSRB.RGCCRP = str2double(fread(fid,7,'uint8=>char')');
MENSRB.RLMAP = fread(fid,1,'uint8=>char')';
MENSRB.RP_ROW = str2double(fread(fid,5,'uint8=>char')');
MENSRB.RP_COL = str2double(fread(fid,5,'uint8=>char')');
MENSRB.C_R_NC = str2double(fread(fid,10,'uint8=>char')');
MENSRB.C_R_EC = str2double(fread(fid,10,'uint8=>char')');
MENSRB.C_R_DC = str2double(fread(fid,10,'uint8=>char')');
MENSRB.C_AZ_NC = str2double(fread(fid,9,'uint8=>char')');
MENSRB.C_AZ_EC = str2double(fread(fid,9,'uint8=>char')');
MENSRB.C_AZ_DC = str2double(fread(fid,9,'uint8=>char')');
MENSRB.C_AL_NC = str2double(fread(fid,9,'uint8=>char')');
MENSRB.C_AL_EC = str2double(fread(fid,9,'uint8=>char')');
MENSRB.C_AL_DC = str2double(fread(fid,9,'uint8=>char')');
MENSRB.TOTAL_TILES_COLS = str2double(fread(fid,3,'uint8=>char')');
MENSRB.TOTAL_TILES_ROWS = str2double(fread(fid,5,'uint8=>char')');

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////