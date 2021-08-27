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
cphd_preamble.Classification = cphdxmeta.CollectionID.Classification;
cphd_preamble.DateTime = cphdxmeta.Global.Timeline.CollectionStart;
cphd_preamble.DataSetID = cphdxmeta.CollectionID.CoreName;
switch upper(cphdxmeta.CollectionID.RadarMode.ModeType)
    case 'SPOTLIGHT'
        cphd_preamble.Mode = 'Spotlight';
    case 'DYNAMIC STRIPMAP'
        cphd_preamble.Mode = 'Scan';
    case 'STRIPMAP'
        cphd_preamble.Mode = 'Strip';
end
if isfield(cphdxmeta.CollectionID,'CollectType')
    cphd_preamble.Geometry = cphdxmeta.CollectionID.CollectType;
end
cphd_preamble.FixedSRP = cphdxmeta.Channel.SRPFixedCPHD;
cphd_preamble.Datum = 'WGS84'; % Only supported value
switch upper(cphdxmeta.Data.SignalArrayFormat)
    case 'CI2'
        cphd_preamble.PHDataType = 'cmplxb';
    case 'CI4'
        cphd_preamble.PHDataType = 'cmplxs';
    case 'CF8'
        cphd_preamble.PHDataType = 'cmplxf';
    otherwise
        % Unrecognized data type.  Default to most precision
        cphd_preamble.PHDataType = 'cmplxf';
end
cphd_preamble.Interleaved = false;
cphd_preamble.TOASaved = cphdxmeta.Channel.Parameters(channel).TOASaved;
cphd_preamble.PhaseSgn = cphdxmeta.Global.SGN;
cphd_preamble.Sensor = cphdxmeta.CollectionID.CollectorName;
cphd_preamble.Grid = 'Polar'; % All CPHDX is of this type
cphd_preamble.DeskewApplied	= true; % CPHDX assumes that deskew has been applied if appropriate
cphd_preamble.FreqReferenceIndex = isfield(cphdxmeta.Global,'RefFreqIndex') && cphdxmeta.Global.RefFreqIndex;
cphd_preamble.NominalCenterFreq = cphdxmeta.Channel.Parameters(channel).FxC;

if isfield(cphdxmeta,'TxRcv')
    if isfield(cphdxmeta.TxRcv,'TxWFParameters')
        if isfield(cphdxmeta.TxRcv.TxWFParameters,'LFMRate')
            cphd_preamble.NominalChirpRate = cphdxmeta.TxRcv.TxWFParameters.LFMRate;
        else
            cphd_preamble.NominalChirpRate = 0;
        end
        if isfield(cphdxmeta.TxRcv.TxWFParameters,'PulseLength')
            cphd_preamble.XmitPulseDuration = cphdxmeta.TxRcv.TxWFParameters.PulseLength;
        end
    end
    if isfield(cphdxmeta.TxRcv,'RcvParameters') && isfield(cphdxmeta.TxRcv.RcvParameters,'SampleRate')
        cphd_preamble.NominalADRate = cphdxmeta.TxRcv.RcvParameters.SampleRate;
    end
    % IFP4 (which took CPHD version 3.0) at one point required the CPHDs it
    % processed to have the TOASaved value be the line below, even though
    % the actual value should probably be related to swath time (receive
    % window length minus pulse length), which describes the set of
    % time-of-arrivals received, and not the frequency spacing, which
    % describes the width off the FFT'd pulse.
    % cphd_preamble.TOASaved = cphd_preamble.NominalADRate/cphd_preamble.NominalChirpRate;
end

cphd_preamble.Nchannels	= cphdxmeta.Data.NumCPHDChannels;
cphd_preamble.Nvectors = cphdxmeta.Data.Channel(channel).NumVectors;
cphd_preamble.Nsamples = cphdxmeta.Data.Channel(channel).NumSamples;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////