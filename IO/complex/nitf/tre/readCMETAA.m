function  CMETAA  = readCMETAA(fid)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
% READCMETAA Create structure for metadata fields.
%   CMETAA = READCMETAA(FID) returns a structure of fields and their
%   associated values.

% Written by Matt Donath, NGA, matthew.b.donath@nga.ic.gov
% References:
%  > MIL-STD-2500A, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.0
%  > MIL-STD-2500C, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.1
%  > STDI-0001,     NATIONAL SUPPORT DATA EXTENSIONS (SDE)(VERSION 1.3/CN2)
%                   FOR THE NATIONAL IMAGERY TRANSMISSION FORMAT (NITF)
%  > STDI-0002,     THE COMPENDIUM OF CONTROLLED EXTENSIONS (CE) FOR THE 
%                   NATIONAL IMAGERY TRANSMISSION FORMAT (NITF) VERSION 2.1
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Create the metadata structure for CMETAA.

CMETAA.CETAG = fread(fid,6,'uint8=>char')';
CMETAA.CEL = str2double(fread(fid,5,'uint8=>char')');

CMETAA.RELATED_TRES = strtrim(fread(fid,2,'uint8=>char')');
CMETAA.RELATED_TRES = strtrim(fread(fid,120,'uint8=>char')');
CMETAA.RD_PRC_NO = strtrim(fread(fid,12,'uint8=>char')');
CMETAA.IF_PROCESS = strtrim(fread(fid,4,'uint8=>char')');
CMETAA.RD_CEN_FREQ = strtrim(fread(fid,4,'uint8=>char')');
CMETAA.RD_MODE = strtrim(fread(fid,5,'uint8=>char')');
CMETAA.RD_PATCH_NO = str2double(fread(fid,4,'uint8=>char')');
            
CMETAA.CMPLX_DOMAIN = strtrim(fread(fid,5,'uint8=>char')');
CMETAA.CMPLX_MAG_REMAP_TYPE = strtrim(fread(fid,4,'uint8=>char')');
CMETAA.CMPLX_LIN_SCALE = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CMPLX_AVG_POWER = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CMPLX_LINLOG_TP = str2double(fread(fid,5,'uint8=>char')');
CMETAA.CMPLX_PHASE_QUANT_FLAG = strtrim(fread(fid,3,'uint8=>char')');
CMETAA.CMPLX_PHASE_QUANT_BIT_DEPTH = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_SIZE_1 = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_IC_1 = strtrim(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_SIZE_2 = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_IC_2 = strtrim(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_IC_BPP = str2double(fread(fid,5,'uint8=>char')');
CMETAA.CMPLX_WEIGHT = strtrim(fread(fid,3,'uint8=>char')');           
CMETAA.CMPLX_AZ_SLL = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_RNG_SLL = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_AZ_TAY_NBAR = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_RNG_TAY_NBAR = str2double(fread(fid,2,'uint8=>char')');
CMETAA.CMPLX_WEIGHT_NORM = strtrim(fread(fid,3,'uint8=>char')');
CMETAA.CMPLX_SIGNAL_PLANE = strtrim(fread(fid,1,'uint8=>char')');
            
CMETAA.IF_DC_SF_ROW = str2double(fread(fid,6,'uint8=>char')');
CMETAA.IF_DC_SF_COL = str2double(fread(fid,6,'uint8=>char')');

for lp = 1:4
    CMETAA.IF_PATCH_ROW(lp) = str2double(fread(fid,6,'uint8=>char')');
    CMETAA.IF_PATCH_COL(lp) = str2double(fread(fid,6,'uint8=>char')');
end

CMETAA.IF_DC_IS_ROW = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_DC_IS_COL = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_IMG_ROW_DC = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_IMG_COL_DC = str2double(fread(fid,8,'uint8=>char')');

for lp = 1:4
    CMETAA.IF_TILE_ROW(lp) = str2double(fread(fid,6,'uint8=>char')');
    CMETAA.IF_TILE_COL(lp) = str2double(fread(fid,6,'uint8=>char')');
end

CMETAA.IF_RD = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_RGWLK = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_KEYSTN = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_LINSFT = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_SUBPATCH = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_GEODIST = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_RGFO = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_BEAM_COMP = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_RGRES = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_AZRES = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_RSS = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_AZSS = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_RSR = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_AZSR = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_RFFT_SAMP = str2double(fread(fid,7,'uint8=>char')');
CMETAA.IF_AZFFT_SAMP = str2double(fread(fid,7,'uint8=>char')');
CMETAA.IF_RFFT_TOT = str2double(fread(fid,7,'uint8=>char')');
CMETAA.IF_AZFFT_TOT = str2double(fread(fid,7,'uint8=>char')');
CMETAA.IF_SUBP_ROW = str2double(fread(fid,6,'uint8=>char')');
CMETAA.IF_SUBP_COL = str2double(fread(fid,6,'uint8=>char')');
CMETAA.IF_SUB_RG = str2double(fread(fid,4,'uint8=>char')');
CMETAA.IF_SUB_AZ = str2double(fread(fid,4,'uint8=>char')');
CMETAA.IF_RFFTS = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_AFFTS = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_RANGE_DATA = strtrim(fread(fid,7,'uint8=>char')');
CMETAA.IF_INCPH = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.IF_SR_NAME1 = strtrim(fread(fid,8,'uint8=>char')');
CMETAA.IF_SR_AMOUNT1 = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_SR_NAME2 = strtrim(fread(fid,8,'uint8=>char')');
CMETAA.IF_SR_AMOUNT2 = str2double(fread(fid,8,'uint8=>char')');
CMETAA.IF_SR_NAME3 = strtrim(fread(fid,8,'uint8=>char')');
CMETAA.IF_SR_AMOUNT3 = str2double(fread(fid,8,'uint8=>char')');

for lp = 1:3
    CMETAA.AF_TYPE{lp} = strtrim(fread(fid,5,'uint8=>char')');
end  

CMETAA.POL_TR = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_RE = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_REFERENCE = strtrim(fread(fid,40,'uint8=>char')');
CMETAA.POL = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_REG = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_ISO_1 = str2double(fread(fid,5,'uint8=>char')');
CMETAA.POL_BAL = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_BAL_MAG = str2double(fread(fid,8,'uint8=>char')');
CMETAA.POL_BAL_PHS = str2double(fread(fid,8,'uint8=>char')');
CMETAA.POL_HCOMP = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_HCOMP_BASIS = strtrim(fread(fid,10,'uint8=>char')');

for lp = 1:3
    CMETAA.POL_HCOMP_COEF(lp) = str2double(fread(fid,9,'uint8=>char')');
end

CMETAA.POL_AFCOMP = strtrim(fread(fid,1,'uint8=>char')');
CMETAA.POL_SPARE_A = strtrim(fread(fid,15,'uint8=>char')');
CMETAA.POL_SPARE_N = strtrim(fread(fid,9,'uint8=>char')');
            
CMETAA.T_UTC_YYYYMMMDD = strtrim(fread(fid,9,'uint8=>char')');
CMETAA.T_HHMMSSUTC = strtrim(fread(fid,6,'uint8=>char')');
CMETAA.T_HHMMSSLOCAL = strtrim(fread(fid,6,'uint8=>char')');
            
CMETAA.CG_SRAC = str2double(fread(fid,11,'uint8=>char')');
CMETAA.CG_SLANT_CONFIDENCE = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CG_CROSS = str2double(fread(fid,11,'uint8=>char')');
CMETAA.CG_CROSS_CONFIDENCE = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CG_CAAC = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_CONE_CONFIDENCE = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_GPSAC = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_GPSAC_CONFIDENCE = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_SQUINT = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_GAAC = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CG_GAAC_CONFIDENCE = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_INCIDENT = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CG_SLOPE = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CG_TILT = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_LD = fread(fid,1,'uint8=>char')';
CMETAA.CG_NORTH = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_NORTH_CONFIDENCE = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_EAST = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_RLOS = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_LOS_CONFIDENCE = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_LAYOVER = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_SHADOW = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_OPM = str2double(fread(fid,7,'uint8=>char')');
CMETAA.CG_MODEL = strtrim(fread(fid,5,'uint8=>char')');
            
CMETAA.CG_AMPT_X = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_AMPT_Y = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_AMPT_Z = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_AP_CONF_XY = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_AP_CONF_Z = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_APCEN_X = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_APCEN_Y = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_APCEN_Z = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_APER_CONF_XY = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_APER_CONF_Z = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_FPNUV_X = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_FPNUV_Y = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_FPNUV_Z = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_IDPNUVX = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_IDPNUVY = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_IDPNUVZ = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_SCECN_X = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_SCECN_Y = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_SCECN_Z = str2double(fread(fid,13,'uint8=>char')');
CMETAA.CG_SC_CONF_XY = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_SC_CONF_Z = str2double(fread(fid,6,'uint8=>char')');
CMETAA.CG_SWWD = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_SNVEL_X = str2double(fread(fid,10,'uint8=>char')');
CMETAA.CG_SNVEL_Y = str2double(fread(fid,10,'uint8=>char')');
CMETAA.CG_SNVEL_Z = str2double(fread(fid,10,'uint8=>char')');
CMETAA.CG_SNACC_X = str2double(fread(fid,10,'uint8=>char')');
CMETAA.CG_SNACC_Y = str2double(fread(fid,10,'uint8=>char')');
CMETAA.CG_SNACC_Z = str2double(fread(fid,10,'uint8=>char')');
CMETAA.CG_SNATT_ROLL = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_SNATT_PITCH = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_SNATT_YAW = str2double(fread(fid,8,'uint8=>char')');
CMETAA.CG_GTP_X = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_GTP_Y = str2double(fread(fid,9,'uint8=>char')');
CMETAA.CG_GTP_Z = str2double(fread(fid,9,'uint8=>char')');

CMETAA.CG_MAP_TYPE = strtrim(fread(fid,4,'uint8=>char')');

if strcmp(CMETAA.CG_MAP_TYPE,'GEOD')
    CMETAA.CG_PATCH_LATCEN = str2double(fread(fid,11,'uint8=>char')');
    CMETAA.CG_PATCH_LNGCEN = str2double(fread(fid,12,'uint8=>char')');
    CMETAA.CG_PATCH_LTCORUL = str2double(fread(fid,11,'uint8=>char')');
    CMETAA.CG_PATCH_LGCORUL = str2double(fread(fid,12,'uint8=>char')');
    CMETAA.CG_PATCH_LTCORUR = str2double(fread(fid,11,'uint8=>char')');
    CMETAA.CG_PATCH_LGCORUR = str2double(fread(fid,12,'uint8=>char')');
    CMETAA.CG_PATCH_LTCORLR = str2double(fread(fid,11,'uint8=>char')');
    CMETAA.CG_PATCH_LGCORLR = str2double(fread(fid,12,'uint8=>char')');
    CMETAA.CG_PATCH_LTCORLL = str2double(fread(fid,11,'uint8=>char')');
    CMETAA.CG_PATCH_LNGCOLL = str2double(fread(fid,12,'uint8=>char')');
    CMETAA.CG_PATCH_LAT_CONFIDENCE = str2double(fread(fid,9,'uint8=>char')');
    CMETAA.CG_PATCH_LONG_CONFIDENCE = str2double(fread(fid,9,'uint8=>char')');
elseif strcmp(CMETAA.CG_MAP_TYPE,'MGRS')
    CMETAA.CG_MGRS_CENT = strtrim(fread(fid,23,'uint8=>char')');
    CMETAA.CG_MGRSCORUL = strtrim(fread(fid,23,'uint8=>char')');
    CMETAA.CG_MGRSCORUR = strtrim(fread(fid,23,'uint8=>char')');
    CMETAA.CG_MGRSCORLR = strtrim(fread(fid,23,'uint8=>char')');
    CMETAA.CG_MGRCORLL = strtrim(fread(fid,23,'uint8=>char')');
    CMETAA.CG_MGRS_CONFIDENCE = str2double(fread(fid,7,'uint8=>char')');
    CMETAA.CG_MGRS_PAD = strtrim(fread(fid,11,'uint8=>char')');
elseif strcmp(CMETAA.CG_MAP_TYPE,'NA')
    CMETAA.CG_MAP_TYPE_BLANK = strtrim(fread(fid,133,'uint8=>char')');
end
CMETAA.CG_SPARE_A = strtrim(fread(fid,144,'uint8=>char')');
CMETAA.CA_CALPA = str2double(fread(fid,7,'uint8=>char')');
CMETAA.WF_SRTFR = str2double(fread(fid,14,'uint8=>char')');
CMETAA.WF_ENDFR = str2double(fread(fid,14,'uint8=>char')');
CMETAA.WF_CHRPRT = str2double(fread(fid,10,'uint8=>char')');
CMETAA.WF_WIDTH = str2double(fread(fid,9,'uint8=>char')');
CMETAA.WF_CENFRQ = str2double(fread(fid,13,'uint8=>char')');
CMETAA.WF_BW = str2double(fread(fid,13,'uint8=>char')');
CMETAA.WF_PRF = str2double(fread(fid,7,'uint8=>char')');
CMETAA.WF_PRI = str2double(fread(fid,9,'uint8=>char')');
CMETAA.WF_CDP = str2double(fread(fid,7,'uint8=>char')');
CMETAA.WF_NUMBER_OF_PULSES = str2double(fread(fid,9,'uint8=>char')');

CMETAA.VPH_COND = strtrim(fread(fid,1,'uint8=>char')');

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////