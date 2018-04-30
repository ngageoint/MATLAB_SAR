function  AIMIDA  = readAIMIDA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READAIMIDA Create structure for metadata fields.
%   AIMIDA = READAIMIDA(FID) returns a structure of fields and their
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

%% Create the metadata structure for AIMIDA.

AIMIDA.CETAG = fread(fid,6,'uint8=>char')';
AIMIDA.CEL = str2double(fread(fid,5,'uint8=>char')');

AIMIDA.MISSION_DATE = fread(fid,7,'uint8=>char')';
AIMIDA.MISSION_NO = fread(fid,4,'uint8=>char')';
AIMIDA.FLIGHT_NO = fread(fid,2,'uint8=>char')';
AIMIDA.OP_NUM = fread(fid,3,'uint8=>char')';
AIMIDA.START_SEGMENT = fread(fid,2,'uint8=>char')';
AIMIDA.REPRO_NUM = str2double(fread(fid,2,'uint8=>char')');
AIMIDA.REPLAY = fread(fid,3,'uint8=>char')';
AIMIDA.Reserved1 = fread(fid,1,'uint8=>char')';
AIMIDA.START_COLUMN = str2double(fread(fid,2,'uint8=>char')');
AIMIDA.START_ROW = str2double(fread(fid,5,'uint8=>char')');
AIMIDA.END_SEGMENT = fread(fid,2,'uint8=>char')';
AIMIDA.END_COLUMN = str2double(fread(fid,2,'uint8=>char')');
AIMIDA.END_ROW = str2double(fread(fid,5,'uint8=>char')');
AIMIDA.COUNTRY = fread(fid,2,'uint8=>char')';
AIMIDA.Reserved2 = fread(fid,4,'uint8=>char')';
AIMIDA.LOCATION = fread(fid,11,'uint8=>char')';
AIMIDA.TIME = fread(fid,5,'uint8=>char')';
AIMIDA.CREATION_DATE = fread(fid,7,'uint8=>char')';

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////