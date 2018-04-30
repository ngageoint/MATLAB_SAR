function TRE = readUndefinedTRE(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READUNDEFINEDTRE Capture data for TREs that currently do not a metadata
% parser.
%   TRE = readUndefinedTRE(FID) returns a structure of unparsed metadata.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Written by Matt Donath, NGA, matthew.b.donath@nga.ic.gov

%% Create the metadata structure for undefined TREs.

TRE.TRETAG = fread(fid,6,'uint8=>char')';
TRE.TREL = str2double(fread(fid,5,'uint8=>char')');

TRE.TREDATA = fread(fid,TRE.TREL,'uint8=>char')'; 

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////