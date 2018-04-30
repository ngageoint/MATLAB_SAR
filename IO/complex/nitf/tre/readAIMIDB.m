function  AIMIDB  = readAIMIDB(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READAIMIDB Create structure for metadata fields.
%   AIMIDB = READAIMIDB(FID) returns a structure of fields and their
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

%% Create the metadata structure for AIMIDB.

AIMIDB.CETAG = fread(fid,6,'uint8=>char')';
AIMIDB.CEL = str2double(fread(fid,5,'uint8=>char')');

AIMIDB.ACQUISITION_DATE = fread(fid,14,'uint8=>char')';
AIMIDB.MISSION_NO = fread(fid,4,'uint8=>char')';
AIMIDB.MISSION_IDENTIFICATION = fread(fid,10,'uint8=>char')';
AIMIDB.FLIGHT_NO = fread(fid,2,'uint8=>char')';
AIMIDB.OP_NUM = fread(fid,3,'uint8=>char')';
AIMIDB.CURRENT_SEGMENT = fread(fid,2,'uint8=>char')';
AIMIDB.REPRO_NUM = str2double(fread(fid,2,'uint8=>char')');
AIMIDB.REPLAY = fread(fid,3,'uint8=>char')';
AIMIDB.RESERVED001 = fread(fid,1,'uint8=>char')';
AIMIDB.START_TILE_COLUMN = str2double(fread(fid,3,'uint8=>char')');
AIMIDB.START_TILE_ROW = str2double(fread(fid,5,'uint8=>char')');
AIMIDB.END_SEGMENT = fread(fid,2,'uint8=>char')';
AIMIDB.END_TILE_COLUMN = str2double(fread(fid,3,'uint8=>char')');
AIMIDB.END_TILE_ROW = str2double(fread(fid,5,'uint8=>char')');
AIMIDB.COUNTRY = fread(fid,2,'uint8=>char')';
AIMIDB.RESERVED002 = fread(fid,4,'uint8=>char')';
AIMIDB.LOCATION = fread(fid,11,'uint8=>char')';
AIMIDB.RESERVED003 = fread(fid,13,'uint8=>char')';

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////