function [ output_meta ] = meta2sicd_aimidb( aimidb_struct )
%META2SICD_AIMIDB Converts metadata stored in the AIMIDB NITF TRE into a
% SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if any(~isspace(aimidb_struct.COUNTRY)) % Could be blank
    output_meta.CollectionInfo.CountryCode=aimidb_struct.COUNTRY;
end
output_meta.ImageFormation.SegmentIdentifier=aimidb_struct.CURRENT_SEGMENT;
try
    output_meta.Timeline.CollectStart=datenum(aimidb_struct.ACQUISITION_DATE,'yyyymmddHHMMSS');
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////