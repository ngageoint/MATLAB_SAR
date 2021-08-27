function [sicdmeta] = pfa_sicd_meta(meta, nbdata, ifp_params)
%PFA_SICD_META Creates a SICD metadata structure for an image created by pfa_file.m
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
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Create SICD meatdata structure
sicdmeta = meta2sicd_cphdx(meta,nbdata,ifp_params.channel);

% Image Creation
sicdmeta.ImageCreation.DateTime = now();
sicdmeta.ImageCreation.Profile  = 'pfa_sicd_meta.m';

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
center=mean(nbdata.SRPPos); % Current PFA code only forms to this point
sicdmeta.GeoData.SCP.ECF.X=mean(center(1));
sicdmeta.GeoData.SCP.ECF.Y=mean(center(2));
sicdmeta.GeoData.SCP.ECF.Z=mean(center(3));
scp_lla = ecf_to_geodetic([sicdmeta.GeoData.SCP.ECF.X sicdmeta.GeoData.SCP.ECF.Y sicdmeta.GeoData.SCP.ECF.Z]);
sicdmeta.GeoData.SCP.LLH.Lat=scp_lla(1);
sicdmeta.GeoData.SCP.LLH.Lon=scp_lla(2);
sicdmeta.GeoData.SCP.LLH.HAE=scp_lla(3);
% Corner coordinates will be popupated later in derived_sicd_fields call

% Grid
sicdmeta.Grid.ImagePlane     = 'SLANT';
sicdmeta.Grid.Type           = 'RGAZIM';
sicdmeta.Grid.Row.Sgn=-1; % We have done an IFFT in PFA
sicdmeta.Grid.Row.ImpRespBW  = diff(ifp_params.k_v_bounds);
sicdmeta.Grid.Row.ImpRespWid = 0.886/sicdmeta.Grid.Row.ImpRespBW;
sicdmeta.Grid.Row.SS         = 1/(ifp_params.sample_rate*sicdmeta.Grid.Row.ImpRespBW);
sicdmeta.Grid.Row.KCtr       = mean(ifp_params.k_v_bounds);
sicdmeta.Grid.Row.DeltaKCOAPoly = 0;
sicdmeta.Grid.Row.WgtType.WindowName = 'UNIFORM';
sicdmeta.Grid.Col.Sgn=-1; % We have done an IFFT in PFA
sicdmeta.Grid.Col.ImpRespBW  = diff(ifp_params.k_u_bounds);
sicdmeta.Grid.Col.ImpRespWid = 0.886/sicdmeta.Grid.Col.ImpRespBW;
sicdmeta.Grid.Col.SS         = 1/(ifp_params.sample_rate*sicdmeta.Grid.Col.ImpRespBW);
sicdmeta.Grid.Col.KCtr       = 0;
sicdmeta.Grid.Col.WgtType.WindowName = 'UNIFORM';
sicdmeta.Grid.Col.DeltaKCOAPoly = 0;

% ImageFormation
sicdmeta.ImageFormation.RcvChanProc.NumChanProc = 1; % This function only handles single-channel data
sicdmeta.ImageFormation.RcvChanProc.ChanIndex = ifp_params.channel;
sicdmeta.ImageFormation.ImageFormAlgo = 'PFA';
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
sicdmeta.ImageFormation.ImageFormAlgo = 'PFA';
sicdmeta.ImageFormation.STBeamComp = 'NO';
sicdmeta.ImageFormation.ImageBeamComp = 'NO';
sicdmeta.ImageFormation.AzAutofocus = 'NO';
sicdmeta.ImageFormation.RgAutofocus = 'NO';

% SCPCOA
sicdmeta.SCPCOA.SCPTime = nbdata.TxTime(ifp_params.ref_pulse_index);

% PFA
sicdmeta.PFA.FPN = struct('X',ifp_params.fpn(1),'Y',ifp_params.fpn(2),'Z',ifp_params.fpn(3));
sicdmeta.PFA.IPN = struct('X',ifp_params.ipn(1),'Y',ifp_params.ipn(2),'Z',ifp_params.ipn(3));
sicdmeta.PFA.PolarAngRefTime = nbdata.TxTime(ifp_params.ref_pulse_index);
% This often givens a "badly conditioned" polynomial warning with fitting a
% polynomial in ECF space.  Ignore.
old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
sicdmeta.PFA.PolarAngPoly = fliplr( polyfit(nbdata.TxTime, ifp_params.k_a, 5) ).';
sicdmeta.PFA.SpatialFreqSFPoly = fliplr( polyfit(ifp_params.k_a, ifp_params.k_sf, 5) ).';
warning(old_state);
sicdmeta.PFA.Krg1 = ifp_params.k_v_bounds(1);
sicdmeta.PFA.Krg2 = ifp_params.k_v_bounds(2);
sicdmeta.PFA.Kaz1 = ifp_params.k_u_bounds(1);
sicdmeta.PFA.Kaz2 = ifp_params.k_u_bounds(2);

sicdmeta = derived_sicd_fields(sicdmeta); % Computes many derived fields

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////