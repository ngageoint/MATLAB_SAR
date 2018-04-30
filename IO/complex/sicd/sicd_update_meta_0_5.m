function [ sicd_meta_current ] = sicd_update_meta_0_5( sicd_meta_0_5 )
%SICD_UPDATE_META_0_5 Update a SICD metadata structure from version 0.5 to
%current version (whatever that may be)
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

sicd_meta_1_0 = sicd_meta_0_5;

% Add RadarCollection.TxPolarization, now required, but optional prior to version 1.0
if isfield(sicd_meta_0_5, 'RadarCollection') && ...
        ~isfield(sicd_meta_0_5.RadarCollection, 'TxPolarization') && ...
        isfield(sicd_meta_0_5.RadarCollection, 'RcvChannels') && ...
        isfield(sicd_meta_0_5.RadarCollection.RcvChannels, 'ChanParameters')
    if numel(sicd_meta_0_5.RadarCollection.RcvChannels.ChanParameters)==1
        sicd_meta_1_0.RadarCollection.TxPolarization = sicd_meta_0_5.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization(1);
    else
        sicd_meta_1_0.RadarCollection.TxPolarization = 'SEQUENCE';
        tx_pols = {sicd_meta_0_5.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization};
        tx_pols = unique(cellfun(@(x) x(1), tx_pols));
        for i = 1:numel(tx_pols)
            sicd_meta_1_0.RadarCollection.TxSequence.TxStep(i).TxPolarization = tx_pols(i);
        end
        % Note: If there are multiple waveforms and multiple polarizations,
        % there is no deconfliction done here.
    end
end

% RadarCollection.Area.Corner was optional in version 0.5, but required in
% version 1.0.  Fortunately, Corner is easily derived from Plane.
if isfield(sicd_meta_0_5, 'RadarCollection') && ...
        isfield(sicd_meta_0_5.RadarCollection, 'Area') && ...
        ~isfield(sicd_meta_0_5.RadarCollection.Area, 'Corner') && ...
        isfield(sicd_meta_0_5.RadarCollection.Area, 'Plane')
    try % If Plane substructure is misformed, this may fail
        plane = sicd_meta_0_5.RadarCollection.Area.Plane; % For concise notation
        ref_pt = [plane.RefPt.ECF.X plane.RefPt.ECF.Y plane.RefPt.ECF.Z];
        x_uvect = [plane.XDir.UVectECF.X plane.XDir.UVectECF.Y plane.XDir.UVectECF.Z];
        y_uvect = [plane.YDir.UVectECF.X plane.YDir.UVectECF.Y plane.YDir.UVectECF.Z];
        x_offsets = [plane.XDir.FirstLine plane.XDir.FirstLine ...
            plane.XDir.NumLines plane.XDir.NumLines];
        y_offsets = [plane.YDir.FirstSample plane.YDir.NumSamples ...
            plane.YDir.NumSamples plane.YDir.FirstSample];
        for i = 1:4
            acp = ref_pt + x_uvect * plane.XDir.LineSpacing * double(x_offsets(i) - plane.RefPt.Line) + ...
                y_uvect * plane.YDir.SampleSpacing * double(y_offsets(i) - plane.RefPt.Sample);
            sicd_meta_1_0.RadarCollection.Area.Corner.ACP(i) = ...
                cell2struct(num2cell(ecf_to_geodetic(acp)),{'Lat','Lon','HAE'});
        end
    end
end

% PolarizationHVAnglePoly no longer a valid field in version 1.0.
if isfield(sicd_meta_0_5, 'RadarCollection') && ...
        isfield(sicd_meta_0_5.RadarCollection, 'PolarizationHVAnglePoly')
    sicd_meta_1_0.RadarCollection = rmfield(sicd_meta_1_0.RadarCollection, 'PolarizationHVAnglePoly');
end

% Antenna.Tx/Rcv/TwoWay.HPBW no longer a valid field in version 1.0.
if isfield(sicd_meta_0_5, 'Antenna')
    if isfield(sicd_meta_0_5.Antenna, 'Tx') && ...
            isfield(sicd_meta_0_5.Antenna.Tx, 'HPBW')
        sicd_meta_1_0.Antenna.Tx = rmfield(sicd_meta_0_5.Antenna.Tx, 'HPBW');
    end
    if isfield(sicd_meta_0_5.Antenna, 'Rcv') && ...
            isfield(sicd_meta_0_5.Antenna.Rcv, 'HPBW')
        sicd_meta_1_0.Antenna.Rcv = rmfield(sicd_meta_0_5.Antenna.Rcv, 'HPBW');
    end
    if isfield(sicd_meta_0_5.Antenna, 'TwoWay') && ...
            isfield(sicd_meta_0_5.Antenna.TwoWay, 'HPBW')
        sicd_meta_1_0.Antenna.TwoWay = rmfield(sicd_meta_0_5.Antenna.TwoWay, 'HPBW');
    end
end

% NoiseLevel got its own substructure between SICD 0.5 and SICD 1.0
if isfield(sicd_meta_0_5, 'Radiometric') && ...
        isfield(sicd_meta_0_5.Radiometric, 'NoisePoly')
    sicd_meta_1_0.Radiometric.NoiseLevel.NoisePoly = ...
            sicd_meta_0_5.Radiometric.NoisePoly;
    sicd_meta_1_0.Radiometric = rmfield(sicd_meta_1_0.Radiometric, 'NoisePoly');
    if isfield(sicd_meta_0_5.Radiometric, 'NoiseLevelType')
        sicd_meta_1_0.Radiometric.NoiseLevel.NoiseLevelType = ...
            sicd_meta_0_5.Radiometric.NoiseLevelType;
    else
        % Even if NoiseLevelType wasn't given, we know that relative noise
        % levels should be 1 at SCP.
        if abs(sicd_meta_0_5.Radiometric.NoisePoly(1)-1)<eps
            sicd_meta_1_0.Radiometric.NoiseLevel.NoiseLevelType = 'RELATIVE';
        else
            sicd_meta_1_0.Radiometric.NoiseLevel.NoiseLevelType = 'ABSOLUTE';
        end
    end
end
if isfield(sicd_meta_0_5, 'Radiometric') && ...
        isfield(sicd_meta_0_5.Radiometric, 'NoiseLevelType')
    sicd_meta_1_0.Radiometric = rmfield(sicd_meta_1_0.Radiometric, 'NoiseLevelType');
end

% MatchInfo
if isfield(sicd_meta_0_5, 'MatchInfo')
    sicd_meta_1_0.MatchInfo = struct(); % Clear this out so we can reconstruct it
    if isfield(sicd_meta_0_5.MatchInfo.Collect, 'MatchType') % MatchType was optional field in 0.5
        types = unique({sicd_meta_0_5.MatchInfo.Collect.MatchType});
    else
        types = {''}; % TypeID (equivalent of MatchType) required in 1.0
    end
    sicd_meta_1_0.MatchInfo.NumMatchTypes = numel(types);
    sicd_meta_1_0.MatchInfo.MatchType = struct();
    for i=1:numel(sicd_meta_1_0.MatchInfo.NumMatchTypes)
        if isfield(sicd_meta_0_5.MatchInfo.Collect, 'MatchType')
            ind = find(strcmp({sicd_meta_0_5.MatchInfo.Collect.MatchType},types{i}));
        else
            ind = 1:numel(sicd_meta_0_5.MatchInfo.Collect);
        end
        sicd_meta_1_0.MatchInfo.MatchType(i).TypeID = strtrim(types{i});
        sicd_meta_1_0.MatchInfo.MatchType(i).NumMatchCollections = ...
            numel(ind) - 1; % 0.5 included current instance as one of the collections
        for j = 1:numel(ind)
            if isfield(sicd_meta_0_5.MatchInfo.Collect(ind(j)), 'Parameter')
                if iscell(sicd_meta_0_5.MatchInfo.Collect(ind(j)).Parameter) % Multiple parameters
                    current_index = find(cellfun(@(x) strcmpi(strtrim(x.name),'CURRENT_INSTANCE'), ...
                        sicd_meta_0_5.MatchInfo.Collect(ind(j)).Parameter),1);
                    if ~isempty(current_index)
                        current_index = str2double(sicd_meta_0_5.MatchInfo.Collect(ind(j)).Parameter{current_index}.value);
                    end
                elseif strcmpi(strtrim(sicd_meta_0_5.MatchInfo.Collect(ind(j)).Parameter.name), 'CURRENT_INSTANCE') % Single parameter
                    current_index = str2double(sicd_meta_0_5.MatchInfo.Collect(ind(j)).Parameter.value);
                else
                    current_index = [];
                end
            end
            if ~isempty(current_index)
                sicd_meta_1_0.MatchInfo.MatchType(i).CurrentIndex = current_index;
            else
                sicd_meta_1_0.MatchInfo.MatchType(i).MatchCollection(j).CoreName = ...
                    sicd_meta_0_5.MatchInfo.Collect(ind(j)).CoreName;
                if isfield(sicd_meta_0_5.MatchInfo.Collect(ind(j)), 'Parameter')
                    sicd_meta_1_0.MatchInfo.MatchType(i).MatchCollection(j).Parameter = ...
                        sicd_meta_0_5.MatchInfo.Collect(ind(j)).Parameter;
                end
            end
        end
    end
end

% Add AzimAng and LayoverAng to SCPCOA
sicd_meta_1_0 = derived_sicd_fields(sicd_meta_1_0);

% We are now updated to version 1.0.  Now do rest of updates.
% sicd_meta_current = sicd_update_meta_1_0(sicd_meta_1_0); % If we needed
% to to more updates past version 1.0, they would go here
sicd_meta_current = sicd_meta_1_0; % No further updates

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////