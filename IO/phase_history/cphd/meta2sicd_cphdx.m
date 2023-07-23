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

sicdmeta.CollectionInfo = cphdmeta.CollectionID;

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
warning(old_state);

% RadarCollection
sicdmeta.RadarCollection.TxFrequency.Min = mean(nbdata.FX1);
sicdmeta.RadarCollection.TxFrequency.Max = mean(nbdata.FX2);
if isfield(cphdmeta.Channel.Parameters,'Polarization')
    sicdmeta.RadarCollection.TxPolarization = ...
        cphdmeta.Channel.Parameters(channel).Polarization.TxPol;
    sicdmeta.RadarCollection.RcvChannels.ChanParameters.TxRcvPolarization = ...
        [cphdmeta.Channel.Parameters(channel).Polarization.TxPol ':' cphdmeta.Channel.Parameters(channel).Polarization.RcvPol];
end


end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////