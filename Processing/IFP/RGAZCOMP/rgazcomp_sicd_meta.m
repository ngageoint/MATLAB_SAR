function [sicdmeta] = rgazcomp_sicd_meta(meta, nbdata, ifp_params)
%RGAZCOMP_SICD_META Creates a SICD metadata structure for an image created by rgzcomp_file.m
%
% Inputs:
%     meta:          Collection narrowband data.  This is the same
%                    structure returned by the get_meta function for a
%                    object returned by open_ph_reader.
%     nbdata:        CPHD-style per-pulse narrowband data.  This is the
%                    same structure returned as the second return parameter
%                    of the read_cphd function.
%     ifp_params:    Structure of some additional parameters computed
%                    during image formation.
%
% Outputs:
%     sicdmeta:      Image formation and collection information in a
%                    structure compatible with the SICD writer. 
%
% Authors: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Create SICD meatdata structure
sicdmeta = meta2sicd_cphdx(meta,nbdata,ifp_params.channel);

% Image Creation
sicdmeta.ImageCreation.DateTime = now();
sicdmeta.ImageCreation.Profile  = 'rgazcomp_sicd_meta.m';

% ImageData
sicdmeta.ImageData.PixelType               = 'RE32F_IM32F';
sicdmeta.ImageData.NumRows                 = ifp_params.image_size_pixels(1);
sicdmeta.ImageData.NumCols                 = ifp_params.image_size_pixels(2);
sicdmeta.ImageData.FirstRow                = 0;
sicdmeta.ImageData.FirstCol                = 0;
sicdmeta.ImageData.FullImage.NumRows       = sicdmeta.ImageData.NumRows;
sicdmeta.ImageData.FullImage.NumCols       = sicdmeta.ImageData.NumCols;
sicdmeta.ImageData.SCPPixel.Row            = floor(sicdmeta.ImageData.NumRows/2);
sicdmeta.ImageData.SCPPixel.Col            = floor(sicdmeta.ImageData.NumCols/2);

% GeoData
sicdmeta.GeoData.EarthModel='WGS_84';
center=mean(nbdata.SRPPos); % IFP assumes a constant SRPPos
sicdmeta.GeoData.SCP.ECF.X=mean(center(1));
sicdmeta.GeoData.SCP.ECF.Y=mean(center(2));
sicdmeta.GeoData.SCP.ECF.Z=mean(center(3));
scp_lla = ecf_to_geodetic([sicdmeta.GeoData.SCP.ECF.X sicdmeta.GeoData.SCP.ECF.Y sicdmeta.GeoData.SCP.ECF.Z]);
sicdmeta.GeoData.SCP.LLH.Lat=scp_lla(1);
sicdmeta.GeoData.SCP.LLH.Lon=scp_lla(2);
sicdmeta.GeoData.SCP.LLH.HAE=scp_lla(3);
% Corner coordinates will be popupated later in derived_sicd_fields call

% ImageFormation
sicdmeta.ImageFormation.RcvChanProc.NumChanProc = 1; % This function only handles single-channel data
sicdmeta.ImageFormation.RcvChanProc.ChanIndex = ifp_params.channel;
sicdmeta.ImageFormation.TStartProc = nbdata.TxTime(1);
sicdmeta.ImageFormation.TEndProc = nbdata.TxTime(end);
if isfield(meta.Channel,'Parameters') && ...
        isfield(meta.Channel.Parameters(ifp_params.channel),'Polarization') && ...
        all(isfield(meta.Channel.Parameters(ifp_params.channel).Polarization,{'TxPol','RcvPol'}))
    sicdmeta.ImageFormation.TxRcvPolarizationProc = ...
        [meta.Channel.Parameters(ifp_params.channel).Polarization.TxPol ':' ...
            meta.Channel.Parameters(ifp_params.channel).Polarization.RcvPol];
end
sicdmeta.ImageFormation.TxFrequencyProc.MinProc = ...
    min(nbdata.SC0 + (nbdata.SCSS * double(min(ifp_params.sample_range)-1)));
sicdmeta.ImageFormation.TxFrequencyProc.MaxProc = ...
    max(nbdata.SC0 + (nbdata.SCSS * double(max(ifp_params.sample_range)-1)));
sicdmeta.ImageFormation.ImageFormAlgo = 'RGAZCOMP';
sicdmeta.ImageFormation.STBeamComp = 'NO';
sicdmeta.ImageFormation.ImageBeamComp = 'NO';
sicdmeta.ImageFormation.AzAutofocus = 'NO';
sicdmeta.ImageFormation.RgAutofocus = 'NO';

% Grid
sicdmeta.Grid.ImagePlane     = 'SLANT';
sicdmeta.Grid.Type           = 'RGAZIM';
v_coa = (1 + length(nbdata.TxTime))/2;
[cV2T,s,mu] = polyfit((1:length(nbdata.TxTime)).',nbdata.TxTime,5);
sicdmeta.Grid.TimeCOAPoly = polyval(cV2T,v_coa,[],mu);

% We assume constant Fx0 and Fx_SS across all vectors here (checked in rgazcomp_file.m)
s_coa = mean(double([min(ifp_params.sample_range), max(ifp_params.sample_range)]));
fx_coa = nbdata.SC0(round(v_coa)) + nbdata.SCSS(round(v_coa)) * (s_coa - 1);
krg_coa = (2/SPEED_OF_LIGHT) * fx_coa;
% Assumes sample_range is regularly spaced samples (checked in rgazcomp_file.m)
fx_ss_sf = mean(diff(ifp_params.sample_range)); % Scale factor for sample spacing and Fx_SS
krg_ss = (2/SPEED_OF_LIGHT) * nbdata.SCSS(round(v_coa)) * fx_ss_sf;

% Grid.Row
sicdmeta.Grid.Row.Sgn=-1; % We use an IFFT
sicdmeta.Grid.Row.ImpRespBW = double(max(ifp_params.sample_range)-min(ifp_params.sample_range))*krg_ss; 
sicdmeta.Grid.Row.ImpRespWid = 0.886/sicdmeta.Grid.Row.ImpRespBW;
sicdmeta.Grid.Row.SS = 1/((sicdmeta.ImageData.NumRows-1)*krg_ss);
sicdmeta.Grid.Row.KCtr = krg_coa; % s_0 = s_coa in our case
sicdmeta.Grid.Row.DeltaKCOAPoly = 0;
sicdmeta.Grid.Row.WgtType.WindowName = 'UNIFORM';

% SCPCOA and RGAZCOMP
sicdmeta = derived_sicd_fields(sicdmeta);

% Grid.Col
cT2V = polyfit(nbdata.TxTime-sicdmeta.Grid.TimeCOAPoly,(0:(length(nbdata.TxTime)-1)).',5);
st_rate_coa = cT2V(end-1); % MATLAB polynomials have linear term second from end
kaz_ss = krg_coa * ... % Delta_Kaz per Delta_Kv without the LOOK term
    (norm([sicdmeta.SCPCOA.ARPVel.X sicdmeta.SCPCOA.ARPVel.Y sicdmeta.SCPCOA.ARPVel.Z]) * ...
    sind(sicdmeta.SCPCOA.DopplerConeAng) / sicdmeta.SCPCOA.SlantRange) / ...
    st_rate_coa;
sicdmeta.Grid.Col.Sgn=-1; % We use an IFFT
sicdmeta.Grid.Col.ImpRespBW  = (length(nbdata.TxTime)-1)*kaz_ss; 
sicdmeta.Grid.Col.ImpRespWid = 0.886/sicdmeta.Grid.Col.ImpRespBW;
sicdmeta.Grid.Col.SS = 1/((sicdmeta.ImageData.NumCols-1)*kaz_ss);
sicdmeta.Grid.Col.KCtr       = 0; % v_0 = v_coa in our case
sicdmeta.Grid.Col.DeltaKCOAPoly = 0;
sicdmeta.Grid.Col.WgtType.WindowName = 'UNIFORM';

% Grid.Col.DeltaK1/2 and GeoData.ImageCorners
sicdmeta = derived_sicd_fields(sicdmeta);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////