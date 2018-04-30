function  MPDSRA  = readMPDSRA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READMPDSRA Create structure for metadata fields.
%   MPDSRA = READMPDSRA(FID) returns a structure of fields and their
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

%% Create the metadata structure for MPDSRA.

MPDSRA.CETAG = fread(fid,6,'uint8=>char')';
MPDSRA.CEL = str2double(fread(fid,5,'uint8=>char')');

MPDSRA.BLK_NUM = str2double(fread(fid,2,'uint8=>char')');
MPDSRA.IPR = str2double(fread(fid,2,'uint8=>char')');
MPDSRA.NBLKS_IN_WDG = str2double(fread(fid,2,'uint8=>char')');
MPDSRA.ROWS_IN_BLK = str2double(fread(fid,5,'uint8=>char')');
MPDSRA.COLS_IN_BLK = str2double(fread(fid,5,'uint8=>char')');
MPDSRA.ORP_X = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ORP_Y = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ORP_Z = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ORP_ROW = str2double(fread(fid,5,'uint8=>char')');
MPDSRA.ORP_COLUMN = str2double(fread(fid,5,'uint8=>char')');
MPDSRA.FOC_X = str2double(fread(fid,7,'uint8=>char')');
MPDSRA.FOC_Y = str2double(fread(fid,7,'uint8=>char')');
MPDSRA.FOC_Z = str2double(fread(fid,7,'uint8=>char')');
MPDSRA.ARP_TIME = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.RESERVED001 = fread(fid,14,'uint8=>char')';
MPDSRA.ARP_POS_N = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ARP_POS_E = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ARP_POS_D = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ARP_VEL_N = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ARP_VEL_E = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ARP_VEL_D = str2double(fread(fid,9,'uint8=>char')');
MPDSRA.ARP_ACC_N = str2double(fread(fid,8,'uint8=>char')');
MPDSRA.ARP_ACC_E = str2double(fread(fid,8,'uint8=>char')');
MPDSRA.ARP_ACC_D = str2double(fread(fid,8,'uint8=>char')');
MPDSRA.RESERVED002 = fread(fid,13,'uint8=>char')';

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////