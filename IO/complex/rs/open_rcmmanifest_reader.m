function [ readerobj ] = open_rcmmanifest_reader( filename )
%OPEN_RCMMANIFEST_READER Intiates a reader object for RCM SAFE file
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

readerobj = open_rs_reader(fullfile(fileparts(filename), 'metadata', 'product.xml'));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////