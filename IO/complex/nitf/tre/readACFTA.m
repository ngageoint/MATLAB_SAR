function  ACFTA  = readACFTA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READACFTA Create structure for metadata fields.
%   ACFTA = READACFTA(FID) returns a structure of fields and their
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

%%  Create the metadata structure for ACFTA.

ACFTA.CETAG = fread(fid,6,'uint8=>char')';
ACFTA.CEL = str2double(fread(fid,5,'uint8=>char')');

ACFTA.AC_MSN_ID = fread(fid,10,'uint8=>char')';
ACFTA.SCTYPE = fread(fid,1,'uint8=>char')';
ACFTA.SCNUM = fread(fid,4,'uint8=>char')';
ACFTA.SENSOR_ID = fread(fid,3,'uint8=>char')';
ACFTA.PATCH_TOT = str2double(fread(fid,4,'uint8=>char')');
ACFTA.MTI_TOT = str2double(fread(fid,3,'uint8=>char')');
ACFTA.PDATE = fread(fid,7,'uint8=>char')';
ACFTA.IMHOSTNO = str2double(fread(fid,3,'uint8=>char')');
ACFTA.IMREQID = str2double(fread(fid,5,'uint8=>char')');
ACFTA.SCENE_SOURCE = fread(fid,1,'uint8=>char')';
ACFTA.MPLAN = fread(fid,2,'uint8=>char')';
ACFTA.ENTLOC = fread(fid,21,'uint8=>char')';
ACFTA.ENTELV = str2double(fread(fid,6,'uint8=>char')');
ACFTA.EXITLOC = fread(fid,21,'uint8=>char')';
ACFTA.EXITELV = str2double(fread(fid,6,'uint8=>char')');
ACFTA.TMAP = str2double(fread(fid,7,'uint8=>char')');
ACFTA.RCS = str2double(fread(fid,3,'uint8=>char')');
ACFTA.ROW_SPACING = str2double(fread(fid,7,'uint8=>char')');
ACFTA.COL_SPACING = str2double(fread(fid,7,'uint8=>char')');
ACFTA.SENSERIAL = fread(fid,4,'uint8=>char')';
ACFTA.ABSWVER = fread(fid,7,'uint8=>char')';

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////