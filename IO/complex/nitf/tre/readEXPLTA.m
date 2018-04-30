function  EXPLTA  = readEXPLTA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READEXPLTA Create structure for metadata fields.
%   EXPLTA = READEXPLTA(FID) returns a structure of fields and their
%   associated values.

% Written by Matt Donath, NGA, matthew.b.donath@nga.ic.gov
% References:
%  > MIL-STD-2500A, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.0
%  > MIL-STD-2500C, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.1
%  > STDI-0001,     NATIONAL SUPPORT DATA EXTENSIONS (SDE)(VERSION 1.3/CN2)
%                   FOR THE NATIONAL IMAGERY TRANSMISSION FORMAT (NITF)
%  > STDI-0002,     THE COMPENDIUM OF CONTROLLED EXTENSIONS (CE) FOR THE 
%                   NATIONAL IMAGERY TRANSMISSION FORMAT (NITF) VERSION 2.1
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Create the metadata structure for EXPLTA.

EXPLTA.CETAG = fread(fid,6,'uint8=>char')';
EXPLTA.CEL = str2double(fread(fid,5,'uint8=>char')');

EXPLTA.ANGLE_TO_NORTH = str2double(fread(fid,3,'uint8=>char')');
EXPLTA.SQUINT_ANGLE = str2double(fread(fid,3,'uint8=>char')');
EXPLTA.MODE = fread(fid,3,'uint8=>char')';
EXPLTA.RESERVED001 = fread(fid,16,'uint8=>char')';
EXPLTA.GRAZE_ANG = str2double(fread(fid,2,'uint8=>char')');
EXPLTA.SLOPE_ANG = str2double(fread(fid,2,'uint8=>char')');
EXPLTA.POLAR = fread(fid,2,'uint8=>char')';
EXPLTA.NSAMP = str2double(fread(fid,5,'uint8=>char')');
EXPLTA.RESERVED002 = fread(fid,1,'uint8=>char')';
EXPLTA.SEQ_NUM = fread(fid,1,'uint8=>char')';
EXPLTA.PRIME_ID = fread(fid,12,'uint8=>char')';
EXPLTA.PRIME_BE = fread(fid,15,'uint8=>char')';
EXPLTA.RESERVED003 = fread(fid,1,'uint8=>char')';
EXPLTA.N_SEC = str2double(fread(fid,2,'uint8=>char')');
EXPLTA.IPR = str2double(fread(fid,2,'uint8=>char')');
EXPLTA.RESERVED004 = fread(fid,2,'uint8=>char')';
EXPLTA.RESERVED005 = fread(fid,2,'uint8=>char')';
EXPLTA.RESERVED006 = fread(fid,5,'uint8=>char')';
EXPLTA.RESERVED007 = fread(fid,8,'uint8=>char')';

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////