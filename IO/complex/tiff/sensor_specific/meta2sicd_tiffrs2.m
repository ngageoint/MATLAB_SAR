function [ meta, symmetry ] = meta2sicd_tiffrs2( tiff_filename )
%META2SICD_TIFFRS2 Compile SICD metadata for RS2 data package
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[ meta, symmetry ] = meta2sicd_tiffrs( tiff_filename, 'RS2' );

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////