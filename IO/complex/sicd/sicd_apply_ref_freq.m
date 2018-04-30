function [ sicd_meta ] = sicd_apply_ref_freq( sicd_meta, ref_freq )
%SICD_APPLY_REF_FREQ Apply Ref_Freq offset to all relevant SICD metadata
%
% This should include all of the fields possibly affected by
% RadarCollection.RefFreqIndex.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if isfield(sicd_meta,'RadarCollection')
    if isfield(sicd_meta.RadarCollection,'TxFrequency')
        if isfield(sicd_meta.RadarCollection.TxFrequency,'Min')
            sicd_meta.RadarCollection.TxFrequency.Min = ...
                sicd_meta.RadarCollection.TxFrequency.Min + ref_freq;
        end
        if isfield(sicd_meta.RadarCollection.TxFrequency,'Max')
            sicd_meta.RadarCollection.TxFrequency.Max = ...
                sicd_meta.RadarCollection.TxFrequency.Max + ref_freq;
        end
    end
    if isfield(sicd_meta.RadarCollection,'Waveform') &&...
            isfield(sicd_meta.RadarCollection.Waveform,'WFParameters')
        for i = 1:numel(sicd_meta.RadarCollection.Waveform.WFParameters)
            if isfield(sicd_meta.RadarCollection.Waveform.WFParameters(i),'TxFreqStart')
                sicd_meta.RadarCollection.Waveform.WFParameters(i).TxFreqStart = ...
                    sicd_meta.RadarCollection.Waveform.WFParameters(i).TxFreqStart + ref_freq;
            end
            if isfield(sicd_meta.RadarCollection.Waveform.WFParameters(i),'RcvFreqStart')
                sicd_meta.RadarCollection.Waveform.WFParameters(i).RcvFreqStart = ...
                    sicd_meta.RadarCollection.Waveform.WFParameters(i).RcvFreqStart + ref_freq;
            end
        end
    end
    
end
if isfield(sicd_meta,'ImageFormation')&&isfield(sicd_meta.ImageFormation,'TxFrequencyProc')
    if isfield(sicd_meta.ImageFormation.TxFrequencyProc,'MinProc')
        sicd_meta.ImageFormation.TxFrequencyProc.MinProc = ...
            sicd_meta.ImageFormation.TxFrequencyProc.MinProc + ref_freq;
    end
    if isfield(sicd_meta.ImageFormation.TxFrequencyProc,'MaxProc')
        sicd_meta.ImageFormation.TxFrequencyProc.MaxProc = ...
            sicd_meta.ImageFormation.TxFrequencyProc.MaxProc + ref_freq;
    end
end
if isfield(sicd_meta,'Antenna')
    if isfield(sicd_meta.Antenna,'Tx')&&isfield(sicd_meta.Antenna.Tx,'FreqZero')
        sicd_meta.Antenna.Tx.FreqZero = sicd_meta.Antenna.Tx.FreqZero + ref_freq;
    end
    if isfield(sicd_meta.Antenna,'Rcv')&&isfield(sicd_meta.Antenna.Rcv,'FreqZero')
        sicd_meta.Antenna.Rcv.FreqZero = sicd_meta.Antenna.Rcv.FreqZero + ref_freq;
    end
    if isfield(sicd_meta.Antenna,'TwoWay')&&isfield(sicd_meta.Antenna.TwoWay,'FreqZero')
        sicd_meta.Antenna.TwoWay.FreqZero = sicd_meta.Antenna.TwoWay.FreqZero + ref_freq;
    end
end
if isfield(sicd_meta,'RMA')&&isfield(sicd_meta.RMA,'INCA')&&...
        isfield(sicd_meta.RMA.INCA,'FreqZero')
    sicd_meta.RMA.INCA.FreqZero = sicd_meta.RMA.INCA.FreqZero + ref_freq;
end
sicd_meta.RadarCollection.RefFreqIndex = 0;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////