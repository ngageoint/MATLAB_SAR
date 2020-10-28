function [CollectTime, Bandwidth] = ComputeTimeBandwidth(meta)
%COMPUTETIMEBANDWIDTH Computes collect time and bandwidth given meta

CollectTime = [];
Bandwidth = [];

%CollectTime
if strcmpi(meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    %first check for IPP TStart and TEnd, if not there then use
    %CollectDuration
    if isfield(meta,'ImageFormation') && isfield(meta.ImageFormation,'TStartProc') && ...
       isfield(meta.ImageFormation,'TEndProc') 
        CollectTime = meta.ImageFormation.TEndProc - meta.ImageFormation.TStartProc;
    elseif isfield(meta,'Timeline') && isfield(meta.Timeline,'IPP') && isfield(meta.Timeline.IPP,'Set') && ...
       isfield(meta.Timeline.IPP.Set,'TStart') && isfield(meta.Timeline.IPP.Set,'TEnd')  
        CollectTime = meta.Timeline.IPP.Set.TEnd - meta.Timeline.IPP.Set.TStart;
    elseif isfield(meta,'Timeline') && isfield(meta.Timeline,'CollectDuration')
        CollectTime = meta.Timeline.CollectDuration;
    end
else
    %determine effective time from resolution and geometry
    SCP = [meta.GeoData.SCP.ECF.X; meta.GeoData.SCP.ECF.Y; meta.GeoData.SCP.ECF.Z];
        
    ARP = [meta.SCPCOA.ARPPos.X; meta.SCPCOA.ARPPos.Y; meta.SCPCOA.ARPPos.Z];
    ARV = [meta.SCPCOA.ARPVel.X; meta.SCPCOA.ARPVel.Y; meta.SCPCOA.ARPVel.Z];
        
    Lambda = 2/meta.Grid.Row.KCtr;
    Theta = meta.Grid.Col.ImpRespBW*Lambda/2;
    R = norm(abs(SCP-ARP));    
    V = norm(ARV);
    rov = (R/V);
    CollectTime = abs((Theta*rov)/sind(meta.SCPCOA.DopplerConeAng));    
end

%Bandwidth
if isfield(meta,'Grid') && isfield(meta.Grid,'Row') && isfield(meta.Grid.Row,'ImpRespBW')
    Bandwidth = meta.Grid.Row.ImpRespBW*SPEED_OF_LIGHT/2;
elseif isfield(meta,'ImageFormation') && isfield(meta.ImageFormation,'TxFrequencyProc') && ...
       isfield(meta.ImageFormation.TxFrequencyProc,'MinProc') && ...
       isfield(meta.ImageFormation.TxFrequencyProc,'MaxProc')
    Bandwidth = meta.ImageFormation.TxFrequencyProc.MaxProc - ...
                meta.ImageFormation.TxFrequencyProc.MinProc;
elseif isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
       isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
       isfield(meta.RadarCollection.Waveform.WFParameters,'TxRFBandwidth')
    Bandwidth = meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth;
end

