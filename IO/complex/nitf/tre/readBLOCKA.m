function  BLOCKA  = readBLOCKA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READBLOCKA Create structure for metadata fields.
%   BLOCKA = READBLOCKA(FID) returns a structure of fields and their
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

%% Create the metadata structure for BLOCKA.

BLOCKA.CETAG = fread(fid,6,'uint8=>char')';
BLOCKA.CEL = str2double(fread(fid,5,'uint8=>char')');

BLOCKA.BLOCK_INSTANCE = str2double(fread(fid,2,'uint8=>char')');
BLOCKA.N_GRAY = str2double(fread(fid,5,'uint8=>char')');
BLOCKA.L_LINES = str2double(fread(fid,5,'uint8=>char')');
BLOCKA.LAYOVER_ANGLE = str2double(fread(fid,3,'uint8=>char')');
BLOCKA.SHADOW_ANGLE = str2double(fread(fid,3,'uint8=>char')');
BLOCKA.RESERVED001 = fread(fid,16,'uint8=>char')';
BLOCKA.FRLC_LOC = fread(fid,21,'uint8=>char')';
BLOCKA.LRLC_LOC = fread(fid,21,'uint8=>char')';
BLOCKA.LRFC_LOC = fread(fid,21,'uint8=>char')';
BLOCKA.FRFC_LOC = fread(fid,21,'uint8=>char')';
BLOCKA.RESERVED002 = fread(fid,5,'uint8=>char')';

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////