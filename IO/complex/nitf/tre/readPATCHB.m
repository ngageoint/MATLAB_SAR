function  PATCHB  = readPATCHB(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READPATCHB Create structure for metadata fields.
%   PATCHB = READPATCHB(FID) returns a structure of fields and their
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

%% Create the metadata structure for PATCHB.

PATCHB.CETAG = fread(fid,6,'uint8=>char')';
PATCHB.CEL = str2double(fread(fid,5,'uint8=>char')');

PATCHB.PAT_NO = str2double(fread(fid,4,'uint8=>char')');
PATCHB.LAST_PAT_FLAG = fread(fid,1,'uint8=>char')';
PATCHB.LNSTRT = str2double(fread(fid,7,'uint8=>char')');
PATCHB.LNSTOP = str2double(fread(fid,7,'uint8=>char')');
PATCHB.AZL = str2double(fread(fid,5,'uint8=>char')');
PATCHB.NVL = str2double(fread(fid,5,'uint8=>char')');
PATCHB.FVL = str2double(fread(fid,3,'uint8=>char')');
PATCHB.NPIXEL = str2double(fread(fid,5,'uint8=>char')');
PATCHB.FVPIX = str2double(fread(fid,5,'uint8=>char')');
PATCHB.FRAME = str2double(fread(fid,3,'uint8=>char')');
PATCHB.UTC = str2double(fread(fid,8,'uint8=>char')');
PATCHB.SHEAD = str2double(fread(fid,7,'uint8=>char')');
PATCHB.GRAVITY = str2double(fread(fid,7,'uint8=>char')');
PATCHB.INS_V_NC = str2double(fread(fid,5,'uint8=>char')');
PATCHB.INS_V_EC = str2double(fread(fid,5,'uint8=>char')');
PATCHB.INS_V_DC = str2double(fread(fid,5,'uint8=>char')');
PATCHB.OFFLAT = str2double(fread(fid,8,'uint8=>char')');
PATCHB.OFFLONG = str2double(fread(fid,8,'uint8=>char')');
PATCHB.TRACK = str2double(fread(fid,3,'uint8=>char')');
PATCHB.GSWEEP = str2double(fread(fid,6,'uint8=>char')');
PATCHB.SHEAR = str2double(fread(fid,8,'uint8=>char')');
PATCHB.BATCH_NO = str2double(fread(fid,6,'uint8=>char')');   

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////