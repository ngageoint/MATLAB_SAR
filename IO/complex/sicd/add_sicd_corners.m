function [ output_meta ] = add_sicd_corners( input_meta )
%ADD_SICD_CORNERS Add corner coords to SICD metadata if they can be
%computed.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

output_meta=input_meta;
if isfield(input_meta,'ImageData')&&...
        all(isfield(input_meta.ImageData,{'NumRows','NumCols'}))
    corner_linesample = double([[1 ; 1] ... % Corners of full image
        [1 ; input_meta.ImageData.NumCols] ...
        [input_meta.ImageData.NumRows ; input_meta.ImageData.NumCols] ...
        [input_meta.ImageData.NumRows ; 1]]);
    corner_latlons=point_slant_to_ground(corner_linesample, input_meta);
    if ~isempty(corner_latlons)
        output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=corner_latlons(1,1);
        output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=corner_latlons(2,1);
        output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=corner_latlons(1,2);
        output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=corner_latlons(2,2);
        output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=corner_latlons(1,3);
        output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=corner_latlons(2,3);
        output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=corner_latlons(1,4);
        output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=corner_latlons(2,4);
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////