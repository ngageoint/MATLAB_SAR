function  PATCHA  = readPATCHA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READPATCHA Create structure for metadata fields.
%   PATCHA = READPATCHA(FID) returns a structure of fields and their
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

%% Create the metadata structure for PATCHA.

PATCHA.CETAG = fread(fid,6,'uint8=>char')';
PATCHA.CEL = str2double(fread(fid,5,'uint8=>char')');

PATCHA.PAT_NO = str2double(fread(fid,4,'uint8=>char')');
PATCHA.LAST_PAT_FLAG = fread(fid,1,'uint8=>char')';
PATCHA.LNSTRT = str2double(fread(fid,7,'uint8=>char')');
PATCHA.LNSTOP = str2double(fread(fid,7,'uint8=>char')');
PATCHA.AZL = str2double(fread(fid,5,'uint8=>char')');
PATCHA.NVL = str2double(fread(fid,5,'uint8=>char')');
PATCHA.FVL = str2double(fread(fid,3,'uint8=>char')');
PATCHA.NPIXEL = str2double(fread(fid,5,'uint8=>char')');
PATCHA.FVPIX = str2double(fread(fid,5,'uint8=>char')');
PATCHA.FRAME = str2double(fread(fid,3,'uint8=>char')');
PATCHA.UTC = str2double(fread(fid,8,'uint8=>char')');
PATCHA.SHEAD = str2double(fread(fid,7,'uint8=>char')');
PATCHA.GRAVITY = str2double(fread(fid,7,'uint8=>char')');
PATCHA.INS_V_NC = str2double(fread(fid,5,'uint8=>char')');
PATCHA.INS_V_EC = str2double(fread(fid,5,'uint8=>char')');
PATCHA.INS_V_DC = str2double(fread(fid,5,'uint8=>char')');
PATCHA.OFFLAT = str2double(fread(fid,8,'uint8=>char')');
PATCHA.OFFLONG = str2double(fread(fid,8,'uint8=>char')');
PATCHA.TRACK = str2double(fread(fid,3,'uint8=>char')');
PATCHA.GSWEEP = str2double(fread(fid,6,'uint8=>char')');
PATCHA.SHEAR = str2double(fread(fid,8,'uint8=>char')');

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////