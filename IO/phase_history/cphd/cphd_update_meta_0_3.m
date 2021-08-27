function [ cphd_meta_current ] = cphd_update_meta_0_3( cphd_meta_0_3 )
%CPHD_UPDATE_META_0_3 Update a CPHD metadata structure from version 0.3 to
%current version (whatever that may be)
%
% Author: Wade Schwartzkopf, NGA/Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

cphd_meta_current = cphd_meta_0_3;

%% Sample type
cphd_meta_current.Data = rmfield(cphd_meta_current.Data,'SampleType');
switch cphd_meta_0_3.Data.SampleType
    case 'RE08I_IM08I'
        cphd_meta_current.Data.SignalArrayFormat='CI2';
        datatype_bytes=1;
    case 'RE16I_IM16I'
        cphd_meta_current.Data.SignalArrayFormat='CI4';
        datatype_bytes=2;
    case 'RE32F_IM32F'
        cphd_meta_current.Data.SignalArrayFormat='CF8';
        datatype_bytes=4;
    otherwise
        error('OPEN_CPHD_READER:UNRECOGNIZED_DATATYPE','Unrecognized data type.');
end

% Calculate per channel parameters
for i = 1:cphd_meta_0_3.Data.NumCPHDChannels
    cphd_meta_current.Data = rmfield(cphd_meta_current.Data,'ArraySize');
    cphd_meta_current.Data.Channel(i) = cphd_meta_0_3.Data.ArraySize(i);
    % CPHD 0.3 just assumed arrays were all in order
    if i==1
        cphd_meta_current.Data.Channel(i).SignalArrayByteOffset = 0;
        cphd_meta_current.Data.Channel(i).PVPArrayByteOffset = 0;
    else
        cphd_meta_current.Data.Channel(i).SignalArrayByteOffset = ...
            cphd_meta_current.Data.Channel(i-1).SignalArrayByteOffset + ...
            (2*datatype_bytes*xml_meta.Data.Channel(i-1).NumSamples*xml_meta.Data.Channel(i-1).NumVectors);
        cphd_meta_current.Data.Channel(i).PVPArrayByteOffset = ...
            cphd_meta_current.Data.Channel(i-1).PVPArrayByteOffset + ...
            (xml_meta.Data.NumBytesPVP*xml_meta.Data.Channel(i).NumVectors);
    end
end

cphd_meta_current = rmfield(cphd_meta_current,'SRP');
cphd_meta_current.Channel.SRPFixedCPHD = strcmpi(cphd_meta_0_3.SRP.SRPType, 'FIXEDPT');
cphd_meta_current = rmfield(cphd_meta_current,'CollectionInfo');
cphd_meta_current.CollectionID = cphd_meta_0_3.CollectionInfo;

%% Adjust frequencies in metadata to be true, not offset values, if
% reference frequency is available.
if isfield(cphd_meta_0_3.Global,'RefFreqIndex')&&cphd_meta_0_3.Global.RefFreqIndex&&exist('cphd_ref_freq','file')
    % Get reference frequency from function the user can place in MATLAB path
    ref_freq=cphd_ref_freq();
    
    % Adjust values in XML data
    for i=1:numel(cphd_meta_0_3.Channel.Parameters)
        cphd_meta_0_3.Channel.Parameters(i).FxCtrNom = cphd_meta_0_3.Channel.Parameters(i).FxCtrNom + ref_freq;
    end
    if isfield(cphd_meta_0_3,'Antenna')&&all(isfield(cphd_meta_0_3.Antenna,{'NumTWAnt','TwoWay'}))
        for i = 1:cphd_meta_0_3.Antenna.NumTWAnt
            cphd_meta_0_3.Antenna.TwoWay(i).FreqZero = cphd_meta_0_3.Antenna.TwoWay(i).FreqZero + ref_freq;
        end
    end
end

%% PVP metadata
% Flatten structure described vector-based parameters to just a single level of field names
cphd_meta_current = rmfield(cphd_meta_current,'VectorParameters');
if strcmp(cphd_meta_0_3.Global.DomainType, 'FX')
    cphd_meta_0_3.VectorParameters.SC0 = cphd_meta_0_3.VectorParameters.FxParameters.Fx0;
    cphd_meta_0_3.VectorParameters.SCSS = cphd_meta_0_3.VectorParameters.FxParameters.Fx_SS;
    cphd_meta_0_3.VectorParameters.FX1 = cphd_meta_0_3.VectorParameters.FxParameters.Fx1;
    cphd_meta_0_3.VectorParameters.FX2 = cphd_meta_0_3.VectorParameters.FxParameters.Fx2;
elseif strcmp(cphd_meta_0_3.Global.DomainType, 'TOA')
    cphd_meta_0_3.VectorParameters.SC0 = cphd_meta_0_3.VectorParameters.FxParameters.DeltaTOA0;
    cphd_meta_0_3.VectorParameters.SCSS = cphd_meta_0_3.VectorParameters.FxParameters.TOA_SS;
else
    error('OPEN_CPHD_READER:UNRECOGNIZED_DOMAIN_TYPE','Unrecognized domain type.');
end
vectorParametersCell = fieldnames(cphd_meta_0_3.VectorParameters);
offset = 0;
for i = 1:numel(vectorParametersCell)
    if isnumeric(cphd_meta_0_3.VectorParameters.(vectorParametersCell{i}))
        cphd_meta_current.PVP.(vectorParametersCell{i})=struct('Offset',offset,...
            'Size',cphd_meta_0_3.VectorParameters.(vectorParametersCell{i})/8,...
            'Format','F8');
        offset = offset + cphd_meta_current.PVP.(vectorParametersCell{i}).Size;
    end
end
cphd_meta_current.Data = rmfield(cphd_meta_current.Data,'NumBytesVBP');
cphd_meta_current.Data.NumBytesPVP = cphd_meta_0_3.Data.NumBytesVBP;

% This brings us up to 1.0.  If newer versions are introduced the updates
% to the next version would be called next.

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////