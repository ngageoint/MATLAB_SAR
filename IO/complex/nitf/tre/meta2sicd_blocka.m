function [ output_meta ] = meta2sicd_blocka( blocka_struct )
%META2SICD_BLOCKA Converts metadata stored in the BLOCKA NITF TREs
% into a SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

corner_names={'FRFC','FRLC','LRLC','LRFC'};

% These will probably not end up in the right order, since we don't know
% the orientation of the image.
for i=1:4
    if isletter(blocka_struct.([corner_names{i} '_LOC'])(1))
        output_meta.GeoData.ImageCorners.ICP.(corner_names{i}).Lat=...
            latlonnum(blocka_struct.([corner_names{i} '_LOC'])(1:10));
        output_meta.GeoData.ImageCorners.ICP.(corner_names{i}).Lon=...
            latlonnum(blocka_struct.([corner_names{i} '_LOC'])(11:21));
    else
        output_meta.GeoData.ImageCorners.ICP.(corner_names{i}).Lat=...
            str2double(blocka_struct.([corner_names{i} '_LOC'])(1:10));
        output_meta.GeoData.ImageCorners.ICP.(corner_names{i}).Lon=...
            str2double(blocka_struct.([corner_names{i} '_LOC'])(11:21));
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////