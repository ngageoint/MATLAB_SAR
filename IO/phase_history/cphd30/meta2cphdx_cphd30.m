function cphdxmeta = meta2cphdx_cphd30(cphd_preamble, nbdata)
%META2CPHD_CPHD30 Converts metadata from read_cphd_preamble into CPHD "X"
%format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% CollectionInfo
cphdxmeta.CollectionInfo.CollectorName = cphd_preamble.Sensor;
if isfield(cphd_preamble,'ID') % For older versions of CPHD (< 3.0)
    cphd_preamble.DataSetID = strtok(cphd_preamble.ID);
end
cphdxmeta.CollectionInfo.CoreName = cphd_preamble.DataSetID;
cphdxmeta.CollectionInfo.CollectType = upper(cphd_preamble.Geometry);
switch upper(cphd_preamble.Mode)
    case 'SPOTLIGHT'
        cphdxmeta.CollectionInfo.RadarMode.ModeType = 'SPOTLIGHT';
    case 'SCAN'
        cphdxmeta.CollectionInfo.RadarMode.ModeType = 'DYNAMIC STRIPMAP';
    case 'STRIP'
        cphdxmeta.CollectionInfo.RadarMode.ModeType = 'STRIPMAP';
end
cphdxmeta.CollectionInfo.Classification = cphd_preamble.Classification;

% Data
switch lower(cphd_preamble.PHDataType)
    case 'cmplxn'
        cphdxmeta.Data.SampleType = 'RE08I_IM08I'; % 4-bit type not supported in CPHDx
    case 'cmplxb'
        cphdxmeta.Data.SampleType = 'RE08I_IM08I';
    case 'cmplxs'
        cphdxmeta.Data.SampleType = 'RE16I_IM16I';
    case 'cmplxf'
        cphdxmeta.Data.SampleType = 'RE32F_IM32F';
    otherwise
        % Unrecognized data type
end
cphdxmeta.Data.NumCPHDChannels = cphd_preamble.Nchannels;
cphdxmeta.Data.NumBytesVBP = sum(cell2mat(struct2cell(cphd_preamble.VectorParameters)));
for i=1:cphdxmeta.Data.NumCPHDChannels
    cphdxmeta.Data.ArraySize(i).NumVectors = cphd_preamble.Nvectors;
    cphdxmeta.Data.ArraySize(i).NumSamples = cphd_preamble.Nsamples;
end

% Global
cphdxmeta.Global.DomainType = 'FX';
if isfield(cphd_preamble,'PhaseSgn')
    cphdxmeta.Global.PhaseSGN = cphd_preamble.PhaseSgn;
else
    cphdxmeta.Global.PhaseSGN = -1; % Assume a default if not given
end
if isfield(cphd_preamble,'FreqReferenceIndex')
    cphdxmeta.Global.RefFreqIndex = cphd_preamble.FreqReferenceIndex;
end
if isfield(cphd_preamble,'DateTime')
    cphdxmeta.Global.CollectStart = cphd_preamble.DateTime;
end
cphdxmeta.Global.CollectDuration = nbdata.TxTime(end) - nbdata.TxTime(1);
cphdxmeta.Global.TxTime1 = nbdata.TxTime(1);
cphdxmeta.Global.TxTime2 = nbdata.TxTime(end);
% ImageArea could go here

% Channel
for i=1:cphdxmeta.Data.NumCPHDChannels
    cphdxmeta.Channel.Parameters(i).SRP_Index = 1; % Only one SRP function in CPHD 3.0
    cphdxmeta.Channel.Parameters(i).NomTOARateSF = 1; % TODO: What is this???
    if isfield(cphd_preamble,'NominalCenterFreq')
        cphdxmeta.Channel.Parameters(i).FxCtrNom =cphd_preamble.NominalCenterFreq;
    else
        cphdxmeta.Channel.Parameters(i).FxCtrNom = ...
            mean(nbdata.Fx0 + (nbdata.FxStepSize*cphd_preamble.Nsamples/2));
    end
    if isfield(cphd_preamble,'XmitPulseDuration')&&isfield(cphd_preamble,'NominalChirpRate')
        cphdxmeta.Channel.Parameters(i).BWSavedNom = ...
            cphd_preamble.XmitPulseDuration*cphd_preamble.NominalChirpRate;
    else
        cphdxmeta.Channel.Parameters(i).BWSavedNom = ...
            mean(nbdata.Fx_SS)*cphd_preamble.Nsamples;
    end
    if isfield(cphd_preamble,'TOASaved')
        cphdxmeta.Channel.Parameters(i).TOASavedNom = cphd_preamble.TOASaved;
    end
end

% SRP
if (isfield(cphd_preamble,'FixedSRP')&&cphd_preamble.FixedSRP)||...
        strcmpi(cphd_preamble.Mode,'SPOTLIGHT')
    cphdxmeta.SRP.SRPType = 'FIXEDPT';
    cphdxmeta.SRP.NumSRPs = 1;
    cphdxmeta.SRP.FIXEDPT.SRPPT.X = mean(nbdata.SRPPos(:,1));
    cphdxmeta.SRP.FIXEDPT.SRPPT.Y = mean(nbdata.SRPPos(:,2));
    cphdxmeta.SRP.FIXEDPT.SRPPT.Z = mean(nbdata.SRPPos(:,3));
else
    cphdxmeta.SRP.SRPType = 'PVTPOLY';
    cphdxmeta.SRP.NumSRPs = 1;
    cphdxmeta.SRP.PVTPOLY.SRPPVTPoly.X = fliplr( polyfit(nbdata.TxTime,nbdata.SRPPos(:,1),4) )';
    cphdxmeta.SRP.PVTPOLY.SRPPVTPoly.Y = fliplr( polyfit(nbdata.TxTime,nbdata.SRPPos(:,2),4) )';
    cphdxmeta.SRP.PVTPOLY.SRPPVTPoly.Z = fliplr( polyfit(nbdata.TxTime,nbdata.SRPPos(:,3),4) )';    
end

% Antenna
% CPHD 3.0 contains no antenna information, so this required structure
% cannot be filled out.

% VectorParameters
cphdxmeta.VectorParameters = cphd_preamble.VectorParameters;
% Update fieldnames to CHPDX standard
% Conversion from CPHD 1.2 to CPHD 3.0
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'Fx1', 'Fx0');
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'RxTime', 'RcvTime');
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'RxPos', 'RcvPos');
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'PulseNumber', 'VectorNumber');
% Conversion from CPHD 3.0 to CPHDX
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'SRP', 'SRPPos');
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'AmpSF0', 'AmpSF');
cphdxmeta.VectorParameters = replace_fieldname(cphdxmeta.VectorParameters, 'FxStepSize', 'Fx_SS');

% This code converts to this point converts to CPHD 0.3.  The following
% update will include updates after that.
cphdxmeta = cphd_update_meta(cphdxmeta);
    
    function structure = replace_fieldname(structure, old_fieldname, new_fieldname)
        if isfield(structure,old_fieldname)
            cell_copy = struct2cell(structure);
            f = fieldnames(structure);
            f{strcmp(old_fieldname,f)} = new_fieldname;
            structure = cell2struct(cell_copy,f);
        end
        % The following version is cleaner code, but changes the order of
        % the fields which is bad for the VectorParameters substructure,
        % in which the order defines how the data is written in the file.
        % if isfield(structure,old_fieldname)
        %     structure.(new_fieldname) = structure.(old_fieldname);
        %     structure = rmfield(structure,old_fieldname);
        % end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////