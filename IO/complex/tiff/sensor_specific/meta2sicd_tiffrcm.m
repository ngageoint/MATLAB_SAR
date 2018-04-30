function [ meta, symmetry ] = meta2sicd_tiffrcm( tiff_filename )
%META2SICD_TIFFRCM Compile SICD metadata for RCM data package
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[ meta, symmetry ] = meta2sicd_tiffrs( tiff_filename, 'RCM' );

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////