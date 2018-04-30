function cphd_preamble = meta2cphd30_cphdx(cphdxmeta, channel)
%META2CPHD30_CPHDX Converts metadata from CPHD "X" XML metadata into CPHD
% version 3.0 preamble format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if (nargin<2)
    channel=1;
end

% Setup CPHD file and write preamble
cphd_preamble.Version = '3.0';
cphd_preamble.Classification = cphdxmeta.CollectionInfo.Classification;
cphd_preamble.DateTime = cphdxmeta.Global.CollectStart;
cphd_preamble.DataSetID = cphdxmeta.CollectionInfo.CoreName;
switch upper(cphdxmeta.CollectionInfo.RadarMode.ModeType)
    case 'SPOTLIGHT'
        cphd_preamble.Mode = 'Spotlight';
    case 'DYNAMIC STRIPMAP'
        cphd_preamble.Mode = 'Scan';
    case 'STRIPMAP'
        cphd_preamble.Mode = 'Strip';
end
if isfield(cphdxmeta.CollectionInfo,'CollectType')
    cphd_preamble.Geometry = cphdxmeta.CollectionInfo.CollectType;
end
cphd_preamble.FixedSRP = strcmp(cphdxmeta.SRP.SRPType, 'FIXEDPT');
cphd_preamble.Datum = 'WGS84'; % Only supported value
cphd_preamble.PHDataType = 'cmplxf';
switch upper(cphdxmeta.Data.SampleType)
    case 'RE08I_IM08I'
        cphd_preamble.PHDataType = 'cmplxb';
    case 'RE16I_IM16I'
        cphd_preamble.PHDataType = 'cmplxs';
    case 'RE32F_IM32F'
        cphd_preamble.PHDataType = 'cmplxf';
    otherwise
        % Unrecognized data type.  Default to most precision
        cphd_preamble.PHDataType = 'cmplxf';
end
cphd_preamble.Interleaved = false;
cphd_preamble.TOASaved = cphdxmeta.Channel.Parameters(channel).TOASavedNom;
cphd_preamble.PhaseSgn = cphdxmeta.Global.PhaseSGN;
cphd_preamble.Sensor = cphdxmeta.CollectionInfo.CollectorName;
cphd_preamble.Grid = 'Polar'; % All CPHDX is of this type
cphd_preamble.DeskewApplied	= true; % CPHDX assumes that deskew has been applied if appropriate
cphd_preamble.FreqReferenceIndex = cphdxmeta.Global.RefFreqIndex;
cphd_preamble.NominalCenterFreq = cphdxmeta.Channel.Parameters(channel).FxCtrNom;

% These values are not supported in CPHD X version 0.3.  However, the
% MATLAB SAR Toolbox framework for handling CPHD X often adds the SICD
% RadarCollection field, since it is impossible to make a valid CPHD 3.0 or
% SICD from the formal spec for CPHD X version 0.3.
if isfield(cphdxmeta,'RadarCollection')&&isfield(cphdxmeta.RadarCollection,'Waveform')&&...
    isfield(cphdxmeta.RadarCollection.Waveform,'WFParameters')
    wfp = cphdxmeta.RadarCollection.Waveform.WFParameters(channel);
    if isfield(wfp,'TxFMRate')
        cphd_preamble.NominalChirpRate = wfp.TxFMRate;
    end
    if isfield(wfp,'ADCSampleRate')
        cphd_preamble.NominalADRate = wfp.ADCSampleRate;
    end
    if isfield(wfp,'TxPulseLength')
        cphd_preamble.XmitPulseDuration = wfp.TxPulseLength;
    end
    % IFP4 (which is almost the only reason we would produce CPHD version
    % 3.0 anymore) requires the CPHDs it processes to have the TOASaved
    % value be this, even though the actual value should probably be
    % related to swath time (receive window length minus pulse length),
    % which describes the set of time-of-arrivals received, and not the
    % frequency spacing, which describes the width off the FFT'd pulse.
    % There is no way to get the TOASaved value that IFP4 wants from CPHDX
    % if the RadarCollection field and required subfields are not defined.
    if isfield(wfp,'TxFMRate')&&isfield(wfp,'ADCSampleRate')
        cphd_preamble.TOASaved = wfp.ADCSampleRate/wfp.TxFMRate;
    end
end

cphd_preamble.Nchannels	= cphdxmeta.Data.NumCPHDChannels;
cphd_preamble.Nvectors = cphdxmeta.Data.ArraySize(channel).NumVectors;
cphd_preamble.Nsamples = cphdxmeta.Data.ArraySize(channel).NumSamples;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////