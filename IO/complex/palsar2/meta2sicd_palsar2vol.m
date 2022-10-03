function [ sicd_meta ] = meta2sicd_palsar2vol( native_meta )
%META2SICD_PALSAR2VOL Converts ALOS PALSAR 2 volume file into a SICD-style metadata structure
%
% Takes as input a metadata structure from read_ceos_vol_meta.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% CollectionInfo
if strncmpi(native_meta.vol_set_id, 'ALOS2', 5)
    sicd_meta.CollectionInfo.CollectorName = 'ALOS2';
    sicd_meta.CollectionInfo.RadarMode.ModeID=native_meta.text.prod_id(9:11);
    % Perhaps ScanSAR modes should be DYNAMIC STRIPMAP?
    if strncmpi('SBS',sicd_meta.CollectionInfo.RadarMode.ModeID,3)
        sicd_meta.CollectionInfo.RadarMode.ModeType='SPOTLIGHT';
    else
        sicd_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
    end
elseif strncmpi(native_meta.vol_set_id, 'STRIX', 5)
    sicd_meta.CollectionInfo.CollectorName = deblank(native_meta.vol_set_id);
    sicd_meta.CollectionInfo.RadarMode.ModeID=native_meta.text.prod_id(9:10);
    if strncmpi('SM',sicd_meta.CollectionInfo.RadarMode.ModeID,2)
        sicd_meta.CollectionInfo.RadarMode.ModeType='STRIPMAP';
    elseif strncmpi('SL',sicd_meta.CollectionInfo.RadarMode.ModeID,2)
        sicd_meta.CollectionInfo.RadarMode.ModeType='DYNAMIC STRIPMAP';
    end
end
% Do we want to convert CoreName to NGA-style pattern?
sicd_meta.CollectionInfo.CoreName = native_meta.text.scene_id(8:end);
sicd_meta.CollectionInfo.CollectType='MONOSTATIC';
sicd_meta.CollectionInfo.Classification='UNCLASSIFIED';

%% ImageCreation
sicd_meta.ImageCreation.DateTime = ...
    datenum([native_meta.log_vol_create_date ...
             native_meta.log_vol_create_time(1:6)] ,'yyyymmddHHMMSS');
sicd_meta.ImageCreation.Profile='Prototype';


end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////