function  EXPLTB  = readEXPLTB(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READEXPLTB Create structure for metadata fields.
%   EXPLTB = READEXPLTB(FID) returns a structure of fields and their
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

%% Create the metadata structure for EXPLTB.

EXPLTB.CETAG = fread(fid,6,'uint8=>char')';
EXPLTB.CEL = str2double(fread(fid,5,'uint8=>char')');

EXPLTB.ANGLE_TO_NORTH = str2double(fread(fid,7,'uint8=>char')');
EXPLTB.ANGLE_TO_NORTH_ACCY = str2double(fread(fid,6,'uint8=>char')');
EXPLTB.SQUINT_ANGLE = str2double(fread(fid,7,'uint8=>char')');
EXPLTB.SQUINT_ANGLE_ACCY = str2double(fread(fid,6,'uint8=>char')');
EXPLTB.MODE = fread(fid,3,'uint8=>char')';
EXPLTB.RESERVED001 = fread(fid,16,'uint8=>char')';
EXPLTB.GRAZE_ANG = str2double(fread(fid,5,'uint8=>char')');
EXPLTB.GRAZE_ANG_ACCY = str2double(fread(fid,5,'uint8=>char')');
EXPLTB.SLOPE_ANG = str2double(fread(fid,5,'uint8=>char')');
EXPLTB.POLAR = fread(fid,2,'uint8=>char')';
EXPLTB.NSAMP = str2double(fread(fid,5,'uint8=>char')');
EXPLTB.RESERVED002 = fread(fid,1,'uint8=>char')';
EXPLTB.SEQ_NUM = fread(fid,1,'uint8=>char')';
EXPLTB.PRIME_ID = fread(fid,12,'uint8=>char')';
EXPLTB.PRIME_BE = fread(fid,15,'uint8=>char')';
EXPLTB.RESERVED003 = fread(fid,1,'uint8=>char')';
EXPLTB.N_SEC = str2double(fread(fid,2,'uint8=>char')');
EXPLTB.IPR = fread(fid,2,'uint8=>char')';
   
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////