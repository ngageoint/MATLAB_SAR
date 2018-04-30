function  ACFTB  = readACFTB(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READACFTB Create structure for metadata fields.
%   ACFTB = READACFTB(FID) returns a structure of fields and their
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

%% Create the metadata structure for ACFTB.

ACFTB.CETAG = fread(fid,6,'uint8=>char')';
ACFTB.CEL = str2double(fread(fid,5,'uint8=>char')');

ACFTB.AC_MSN_ID = fread(fid,20,'uint8=>char')';
ACFTB.AC_TAIL_NO = fread(fid,10,'uint8=>char')';
ACFTB.AC_TO = fread(fid,12,'uint8=>char')';
ACFTB.SENSOR_ID_TYPE = fread(fid,4,'uint8=>char')';
ACFTB.SENSOR_ID = fread(fid,6,'uint8=>char')';
ACFTB.SCENE_SOURCE = fread(fid,1,'uint8=>char')';
ACFTB.SCNUM = str2double(fread(fid,6,'uint8=>char')');
ACFTB.PDATE = fread(fid,8,'uint8=>char')';
ACFTB.IMHOSTNO = fread(fid,6,'uint8=>char')';
ACFTB.IMREQID = fread(fid,5,'uint8=>char')';
ACFTB.MPLAN = fread(fid,3,'uint8=>char')';
ACFTB.ENTLOC = fread(fid,25,'uint8=>char')';
ACFTB.LOC_ACCY = str2double(fread(fid,6,'uint8=>char')');
ACFTB.ENTELV = str2double(fread(fid,6,'uint8=>char')');
ACFTB.ELV_UNIT = fread(fid,1,'uint8=>char')';
ACFTB.EXITLOC = fread(fid,25,'uint8=>char')';
ACFTB.EXITELV = str2double(fread(fid,6,'uint8=>char')');
ACFTB.TMAP = str2double(fread(fid,7,'uint8=>char')');
ACFTB.ROW_SPACING = str2double(fread(fid,7,'uint8=>char')');
ACFTB.ROW_SPACING_UNITS = fread(fid,1,'uint8=>char')';
ACFTB.COL_SPACING = str2double(fread(fid,7,'uint8=>char')');
ACFTB.COL_SPACING_UNITS = fread(fid,1,'uint8=>char')';
ACFTB.FOCAL_LENGTH = str2double(fread(fid,6,'uint8=>char')');
ACFTB.SENSERIAL = fread(fid,6,'uint8=>char')';
ACFTB.ABSWVER = fread(fid,7,'uint8=>char')';
ACFTB.CAL_DATE = fread(fid,8,'uint8=>char')';
ACFTB.PATCH_TOT = str2double(fread(fid,4,'uint8=>char')');
ACFTB.MTI_TOT = str2double(fread(fid,3,'uint8=>char')');

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////