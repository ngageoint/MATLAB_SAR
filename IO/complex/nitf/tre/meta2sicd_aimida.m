function [ output_meta ] = meta2sicd_aimida( aimida_struct )
%META2SICD_AIMIDB Converts metadata stored in the AIMIDA NITF TRE into a
% SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if any(~isspace(aimida_struct.COUNTRY)) % Could be blank
    output_meta.CollectionInfo.CountryCode=aimida_struct.COUNTRY;
end
output_meta.ImageCreation.DateTime=datenum(aimida_struct.CREATION_DATE,'ddmmmyy');
output_meta.Timeline.CollectStart=datenum([aimida_struct.MISSION_DATE aimida_struct.TIME],'ddmmmyyHHMM');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////