function metadata = GetMetaData(filenames,fields)
%GETMETADATA Gets metadata from SAR files for specified fields (when
%available)

metadata = repmat(cell2struct(cell(numel(fields),1),fields),0,0);  % Empty struct

waith = waitbar(0,'Processing Metadata');

for i=1:length(filenames)
    
    Temp = [];
    
    a = sprintf('Processing file %d of %d',i,length(filenames));
    waitbar(double(i)/double(length(filenames)),waith,a);
    PHD = 0;
    try
        reader_obj = open_reader(filenames{i});
        if iscell(reader_obj)
            reader_obj = reader_obj{1};
        end
        meta = reader_obj.get_meta();
        SamplePixel = reader_obj.read_chip([1 1],[1 1]);
        reader_obj.close();
    catch
        try
            reader_obj = open_ph_reader(filenames{i});
            if iscell(reader_obj)
                reader_obj = reader_obj{1};
            end
            meta = reader_obj.get_meta();
            PHD = 1;
            if ~isfield(meta,'SCPCOA')
                % Get info for a few specific pulses
                if isfield(reader_obj,'read_raw')
                    read_fun = reader_obj.read_raw;
                elseif isfield(reader_obj,'read_cphd')
                    read_fun = reader_obj.read_cphd;
                end
                [~, vbmeta] = read_fun([1 ...
                    round(meta.Data.Channel(1).NumVectors/2) ...
                    round(meta.Data.Channel(1).NumVectors/2) + 1 ...
                    meta.Data.Channel(1).NumVectors],[]); %first, last and two center pulses (to get velocity)
            end
            meta = derived_phd_fields(meta,vbmeta);
        catch
            continue;
        end        
    end
    
    %loop through fieldnames and populate as metadat supports
    for j=1:length(fields)
        switch fields{j}
            case 'Corename'
                if isfield(meta,'CollectionInfo') && isfield(meta.CollectionInfo,'CoreName')
                    Temp.Corename = meta.CollectionInfo.CoreName;                    
                elseif isfield(meta,'CollectionID') && isfield(meta.CollectionID,'CoreName')
                    Temp.Corename = meta.CollectionID.CoreName;                    
                else
                    Temp.Corename = '';
                end
            case 'DataType'  %todo: if PHD then set to PHD
                if PHD
                    Temp.DataType = 'PHD';
                else
                    if isreal(SamplePixel)
                        Temp.DataType = 'Detected';
                    else
                        Temp.DataType = 'Complex';
                    end
                end
            case 'PixelType'
                if isfield(meta,'ImageData') && isfield(meta.ImageData,'PixelType')
                    Temp.PixelType = meta.ImageData.PixelType;
                else
                    Temp.PixelType = '';
                end
            case 'Lat'
                if isfield(meta,'GeoData') && isfield(meta.GeoData,'SCP') && ...
                   isfield(meta.GeoData.SCP,'LLH') && isfield(meta.GeoData.SCP.LLH,'Lat')
                    Temp.Lat = meta.GeoData.SCP.LLH.Lat;
                else
                    Temp.Lat = '';
                end
            case 'Lon'
                if isfield(meta,'GeoData') && isfield(meta.GeoData,'SCP') && ...
                   isfield(meta.GeoData.SCP,'LLH') && isfield(meta.GeoData.SCP.LLH,'Lon')
                    Temp.Lon = meta.GeoData.SCP.LLH.Lon;
                else
                    Temp.Lon = '';
                end
            case 'CollectTime'
                if isfield(meta,'Timeline') && isfield(meta.Timeline,'CollectStart')
                    Temp.CollectTime = datestr(meta.Timeline.CollectStart);
                elseif isfield(meta,'Global') && isfield(meta.Global,'CollectStart')
                    Temp.CollectTime = datestr(meta.Global.CollectStart);
                elseif isfield(meta,'Global') && isfield(meta.Global,'Timeline') && ...
                        isfield(meta.Global.Timeline,'CollectionStart')
                    Temp.CollectTime = datestr(meta.Global.Timeline.CollectionStart);
                else
                    Temp.CollectTime = '';
                end
            case 'Azimuth'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'AzimAng')
                    Temp.Azimuth = meta.SCPCOA.AzimAng;
                else
                    Temp.Azimuth = '';
                end
            case 'Graze'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'GrazeAng')
                    Temp.Graze = meta.SCPCOA.GrazeAng;
                else
                    Temp.Graze = '';
                end
            case 'NumRows'
                if isfield(meta,'ImageData') && isfield(meta.ImageData,'NumRows')
                    Temp.NumRows = meta.ImageData.NumRows;
                else
                    Temp.NumRows = '';
                end
            case 'NumColumns'
                if isfield(meta,'ImageData') && isfield(meta.ImageData,'NumCols')
                    Temp.NumColumns = meta.ImageData.NumCols;
                else
                    Temp.NumColumns = '';
                end
            case 'AzSS'
                if isfield(meta,'Grid') && isfield(meta.Grid,'Col') && ...
                   isfield(meta.Grid.Col,'SS')
                    Temp.AzSS = meta.Grid.Col.SS;
                else
                    Temp.AzSS = '';
                end
            case 'RgSS'
                if isfield(meta,'Grid') && isfield(meta.Grid,'Row') && ...
                   isfield(meta.Grid.Row,'SS')
                    Temp.RgSS = meta.Grid.Row.SS;
                else
                    Temp.RgSS = '';
                end    
            case 'Duration' %todo: check for TStart and TEnd
                if isfield(meta,'Timeline') && isfield(meta.Timeline,'CollectDuration')
                    Temp.Duration = meta.Timeline.CollectDuration;
                else
                    Temp.Duration = '';
                end    
            case 'RadarMode'
                if isfield(meta,'CollectionInfo') && isfield(meta.CollectionInfo,'RadarMode') && ...
                   isfield(meta.CollectionInfo.RadarMode,'ModeType')
                    Temp.RadarMode = meta.CollectionInfo.RadarMode.ModeType;
                elseif isfield(meta,'CollectionID') && isfield(meta.CollectionID,'RadarMode') && ...
                   isfield(meta.CollectionID.RadarMode,'ModeType')
                    Temp.RadarMode = meta.CollectionID.RadarMode.ModeType;
                else
                    Temp.RadarMode = '';
                end
            case 'Polarization'
                if isfield(meta,'ImageFormation') && isfield(meta.ImageFormation,'TxRcvPolarizationProc')
                    Temp.Polarization = meta.ImageFormation.TxRcvPolarizationProc;
                elseif isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'RcvChannels') && ...
                       isfield(meta.RadarCollection.RcvChannels,'ChanParameters') && ...
                       isfield(meta.RadarCollection.RcvChannels.ChanParameters(1),'TxRcvPolarization')
                    Temp.Polarization = meta.RadarCollection.RcvChannels.ChanParameters(1).TxRcvPolarization;
                elseif isfield(meta,'Channel') && isfield(meta.Channel,'Parameters') && ...
                        isfield(meta.Channel.Parameters,'Polarization') && ...
                        all(isfield(meta.Channel.Parameters(1).Polarization,{'TxPol','RcvPol'}))
                    Temp.Polarization = [meta.Channel.Parameters(1).Polarization.TxPol ...
                        ':' meta.Channel.Parameters(1).Polarization.RcvPol];
                elseif isfield(meta,'Channel') && isfield(meta.Channel,'Parameters') && ...
                        isfield(meta.Channel.Parameters,'RcvPol')
                    Temp.Polarization = meta.Channel.Parameters(1).RcvPol;
                else                    
                    Temp.Polarization = '';
                end
            case 'ImagePlane'
                if isfield(meta,'Grid') && isfield(meta.Grid,'ImagePlane')
                    Temp.ImagePlane = meta.Grid.ImagePlane;
                else
                    Temp.ImagePlane = '';
                end
            case 'HAE'
                if isfield(meta,'GeoData') && isfield(meta.GeoData,'SCP') && ...
                   isfield(meta.GeoData.SCP,'LLH') && isfield(meta.GeoData.SCP.LLH,'HAE')
                    Temp.HAE = meta.GeoData.SCP.LLH.HAE;
                else
                    Temp.HAE = '';
                end
            case 'SCPTime'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SCPTime')
                    Temp.SCPTime = meta.SCPCOA.SCPTime;
                else
                    Temp.SCPTime = '';
                end      
            case 'Side'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SideOfTrack')
                    Temp.Side = meta.SCPCOA.SideOfTrack;
                else
                    Temp.Side = '';
                end      
            case 'SlantRange'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SlantRange')
                    Temp.SlantRange = meta.SCPCOA.SlantRange;
                else
                    Temp.SlantRange = '';
                end  
            case 'GroundRange'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'GroundRange')
                    Temp.GroundRange = meta.SCPCOA.GroundRange;
                else
                    Temp.GroundRange = '';
                end        
            case 'DCA'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'DopplerConeAng')
                    Temp.DCA = meta.SCPCOA.DopplerConeAng;
                else
                    Temp.DCA = '';
                end   
            case 'Twist'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'TwistAng')
                    Temp.Twist = meta.SCPCOA.TwistAng;
                else
                    Temp.Twist = '';
                end
            case 'Slope'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SlopeAng')
                    Temp.Slope = meta.SCPCOA.SlopeAng;
                else
                    Temp.Slope = '';
                end   
            case 'Layover'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'LayoverAng')
                    Temp.Layover = meta.SCPCOA.LayoverAng;
                else
                    Temp.Layover = '';
                end      
            case 'Filename'
                Temp.Filename = filenames{i};
            case 'ROV'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SlantRange') && ...
                   isfield(meta.SCPCOA,'ARPVel')
                   Vel = norm([meta.SCPCOA.ARPVel.X meta.SCPCOA.ARPVel.Y meta.SCPCOA.ARPVel.Z]); 
                   Temp.ROV = meta.SCPCOA.SlantRange./Vel;
                else
                   Temp.ROV = ''; 
                end    
            case 'PRF'
                if isfield(meta,'Timeline') && isfield(meta.Timeline,'IPP') && isfield(meta.Timeline.IPP,'Set') && ...
                   isfield(meta.Timeline.IPP.Set,'IPPStart') && isfield(meta.Timeline.IPP.Set,'IPPEnd')  && ...
                   isfield(meta.Timeline.IPP.Set,'TStart') && isfield(meta.Timeline.IPP.Set,'TEnd')
                   Temp.PRF = (meta.Timeline.IPP.Set.IPPEnd-meta.Timeline.IPP.Set.IPPStart)./ ...
                                          (meta.Timeline.IPP.Set.TEnd-meta.Timeline.IPP.Set.TStart); 
                else                    
                   Temp.PRF = ''; 
                end
            case 'SARAlt'
                if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'ARPPos')
                    R = norm([meta.SCPCOA.ARPPos.X meta.SCPCOA.ARPPos.Y meta.SCPCOA.ARPPos.Z]);
                    Temp.SARAlt = (R-6378135)./1000; %km 
                elseif isfield(meta,'Position') && isfield(meta.Position,'ARPPoly')
                    R = norm([meta.Position.ARPPoly.X(1) meta.Position.ARPPoly.Y(1) meta.Position.ARPPoly.Z(1)]);
                    Temp.SARAlt = (R-6378135)./1000; %km 
                else
                    Temp.SARAlt = ''; 
                end
            case 'CenterFrequency'
                if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'TxFrequency') && ...
                   all(isfield(meta.RadarCollection.TxFrequency,{'Min','Max'}))
                    Temp.CenterFrequency = (meta.RadarCollection.TxFrequency.Max + ...
                                            meta.RadarCollection.TxFrequency.Min)/2;
                else
                    Temp.CenterFrequency = '';
                end
            case 'TxBandwidth'
                if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'TxFrequency') && ...
                   all(isfield(meta.RadarCollection.TxFrequency,{'Min','Max'}))
                    Temp.TxBandwidth = meta.RadarCollection.TxFrequency.Max - ...
                                            meta.RadarCollection.TxFrequency.Min;
                else
                    Temp.TxBandwidth = '';
                end
            case 'TxPulseLength'
                if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
                   isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
                   isfield(meta.RadarCollection.Waveform.WFParameters,'TxPulseLength')
                    Temp.TxPulseLength = meta.RadarCollection.Waveform.WFParameters(1).TxPulseLength;
                else
                    Temp.TxPulseLength = '';
                end
            case 'TxChirpRate'
                if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
                   isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
                   isfield(meta.RadarCollection.Waveform.WFParameters,'TxFMRate')
                    Temp.TxChirpRate = meta.RadarCollection.Waveform.WFParameters(1).TxFMRate;
                else
                    Temp.TxChirpRate = '';
                end
            case 'SampleRate'
                if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
                   isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
                   isfield(meta.RadarCollection.Waveform.WFParameters,'ADCSampleRate')
                    Temp.SampleRate = meta.RadarCollection.Waveform.WFParameters(1).ADCSampleRate;
                else
                    Temp.SampleRate = '';
                end
            case 'RxWindowLength'
                if isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
                   isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
                   isfield(meta.RadarCollection.Waveform.WFParameters,'RcvWindowLength')
                    Temp.RxWindowLength = meta.RadarCollection.Waveform.WFParameters(1).RcvWindowLength;
                else
                    Temp.RxWindowLength = '';
                end
            case 'NumPulses'
                if isfield(meta,'Data') && isfield(meta.Data,'Channel') && ...
                   isfield(meta.Data.Channel,'NumVectors')
                    Temp.NumPulses = meta.Data.Channel.NumVectors;
                elseif isfield(meta,'Timeline') && isfield(meta.Timeline,'IPP') && ...
                        isfield(meta.Timeline.IPP,'Set') && ...
                        isfield(meta.Timeline.IPP.Set,'IPPEnd')
                    Temp.NumPulses = max([meta.Timeline.IPP.Set.IPPEnd]);
                else
                    Temp.NumPulses = '';
                end
            case 'NumSamples'
                if isfield(meta,'Data') && isfield(meta.Data,'Channel') && ...
                   isfield(meta.Data.Channel,'NumSamples')
                    Temp.NumSamples = meta.Data.Channel.NumSamples;
                elseif isfield(meta,'RadarCollection') && isfield(meta.RadarCollection,'Waveform') && ...
                   isfield(meta.RadarCollection.Waveform,'WFParameters') && ...
                   all(isfield(meta.RadarCollection.Waveform.WFParameters,{'ADCSampleRate','RcvWindowLength'}))
                    Temp.NumSamples = meta.RadarCollection.Waveform.WFParameters(1).ADCSampleRate * ...
                       meta.RadarCollection.Waveform.WFParameters(1).RcvWindowLength;
                else
                    Temp.NumSamples = '';
                end
            otherwise
                Temp.(fields{j}) = '';             
        end        
    end
    metadata(end+1) = Temp;
end

close(waith);


end

