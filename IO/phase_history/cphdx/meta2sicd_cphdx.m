function sicdmeta = meta2sicd_cphdx(cphdmeta, nbdata, channel)
%META2SICD_CPHDX Converts metadata from open_ph_reader into SICD-style structure
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if (nargin<3)
    channel = 1;
end

% Fields that can transfer directly from CPHD to SICD, since they merely
% relate the collection, not the image formation.  Some of these fields may
% not even be official CPHD fields, but we often borrow them from SICD to
% supplement CPHD.
transfer_fname = {'CollectionInfo', ...
    'RadarCollection', ...
    'Antenna', ...
    'ErrorStatistics'};
for i = 1:numel(transfer_fname)
    if isfield(cphdmeta, transfer_fname{i})
        sicdmeta.(transfer_fname{i}) = cphdmeta.(transfer_fname{i});
    end
end

% Timeline
if isfield(cphdmeta.Global,'CollectStart')
    sicdmeta.Timeline.CollectStart = cphdmeta.Global.CollectStart;
end
if isfield(cphdmeta.Global,'CollectDuration')
    sicdmeta.Timeline.CollectDuration = cphdmeta.Global.CollectDuration;
end
% Assumes IPPs can be described by a single polynomial
sicdmeta.Timeline.IPP.Set.TStart = nbdata.TxTime(1);
sicdmeta.Timeline.IPP.Set.TEnd = nbdata.TxTime(end);
sicdmeta.Timeline.IPP.Set.IPPStart = 0;
sicdmeta.Timeline.IPP.Set.IPPEnd = length(nbdata.TxTime)-1;
% Polyfit often givens a "badly conditioned" polynomial warning with
% fitting a polynomial in ECF space.  Ignore in this function.
old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
sicdmeta.Timeline.IPP.Set.IPPPoly = fliplr( ...
    polyfit(nbdata.TxTime,(0:sicdmeta.Timeline.IPP.Set.IPPEnd).',5) ).';

% Position
sicdmeta.Position.ARPPoly.X = fliplr( polyfit(nbdata.TxTime,nbdata.TxPos(:,1),4) )';
sicdmeta.Position.ARPPoly.Y = fliplr( polyfit(nbdata.TxTime,nbdata.TxPos(:,2),4) )';
sicdmeta.Position.ARPPoly.Z = fliplr( polyfit(nbdata.TxTime,nbdata.TxPos(:,3),4) )';
if strcmpi(cphdmeta.SRP.SRPType, 'FIXEDPT')
    sicdmeta.Position.GRPPoly = cphdmeta.SRP.FIXEDPT.SRPPT;
else
    if isfield(nbdata,'SRPTime') % Optional field
        time = nbdata.SRPTime;
    else
        time = nbdata.TxTime; % Required field
    end
    sicdmeta.Position.GRPPoly.X = fliplr( polyfit(time,nbdata.SRPPos(:,1),4) )';
    sicdmeta.Position.GRPPoly.Y = fliplr( polyfit(time,nbdata.SRPPos(:,2),4) )';
    sicdmeta.Position.GRPPoly.Z = fliplr( polyfit(time,nbdata.SRPPos(:,3),4) )';    
end
warning(old_state);

% RadarCollection
% Only populate RadarCollection if we can get real frequency info
if ~isfield(cphdmeta,'RadarCollection') && ...
    ((~isfield(cphdmeta.Global,'RefFreqIndex')) || ...
        (~cphdmeta.Global.RefFreqIndex))
    % This will only be bandwidth processed to CPHD, not necessarily
    % transmitted bandwidth, which these SICD fields were really intended
    % to show.
    sicdmeta.RadarCollection.RefFreqIndex = 0; % Any frequencies we use will be real
    sicdmeta.RadarCollection.TxFrequency.Min = mean(nbdata.Fx0);
    sicdmeta.RadarCollection.Waveform.WFParameters.TxRFBandwidth = ...
        mean(nbdata.Fx_SS)*cphdmeta.Data.ArraySize(channel).NumSamples;
    sicdmeta.RadarCollection.TxFrequency.Max = sicdmeta.RadarCollection.TxFrequency.Min + ...
        sicdmeta.RadarCollection.Waveform.WFParameters.TxRFBandwidth;
    sicdmeta.RadarCollection.Waveform.WFParameters.TxFreqStart = sicdmeta.RadarCollection.TxFrequency.Min;
    % Other waveform fields not supported by CPHDX version 0.3
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////