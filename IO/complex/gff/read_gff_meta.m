function [Header, gffName, gffPath, fid_out] = read_gff_meta( FileName, PathName )
% function <[>Header<, gffName, gffPath <,fid_out>]> = read_gff_meta( <FileName> )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% This software is provided as is, with no guarantee, 
%% warranty, or other assurance of functionality or 
%% correctness for any purpose.  Furthermore, users should
%% have no expectation of support from the authors or 
%% Sandia National Laboratories.  Use at your own risk.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function reads the header of an image file written according 
% to the GSAT File Format specification [1].  It returns the values
% from the header in a MATLAB structure, and (optionally) the file
% handle to the valid GFF file.
%
% It leaves the file open (with the proper "endian" setting), and 
% positioned just past the end of the header (at the beginning of
% pixel data).
%
% If the file is NOT in GFF or is invalid in any other way, it is 
% closed and fid_out will be < 0.  If the requested file does not
% exist, fid_out will be -2.
%
% Inputs:
% FileName (Optional)  If the optional file name argument is present, 
%         it is assumed to be the name of a valid, existing GFF file.
%         
%         If the argument is not present, the user is prompted to
%         select a file using the file open dialog box.
%
% Outputs:
% Header  The MATLAB structure containing all fields of the file header.
% gffName (Optional) The name of the opened file.  The input name is
%         returned if "FileName" is passed in.
% gffPath (Optional) The path to the opened file.  Null string returned
%         if "FileName" is passed in.
% fid_out (Optional) This function will sometimes close the file and
%         re-open it.  It also can be used to find a new file without
%         having a handle passed in.  In either case, the calling 
%         routine can get the active file handle from this parameter. 

% References:
% 1) Mendelsohn, G. H., et al, "GSAT Image File Format Specification"
%    Revision 1-6 (Draft), 02Feb00.

% Author:    W. H. Hensley, 2344
% Written:   4Feb2000
% Copyright: 2000 Sandia National Laboratories

% Revisions:
% 08Feb2000 WHH Changed scaling back to 2^-16 on APCAlt.
% 09Feb2000 WHH Changed scaling on RgChirpRate to 2^+12 from 2^-16!
%               Scaling on RgChirpRate now dependent on rev number.
% 27May2000 WHH Changed input argument to FileName and implemented.
% 08Jul2000 DLB/WHH Added multiple-extension support for Unix
%               Added default path support.
% 18Jul2000 WHH Now closes file if the 4th input parameter (fid_out)
%               is NOT supplied.
% 21Aug2000 WHH/DLB Fixed the order of AzResolution and RgResolution.
%               Prior to this, function did NOT meet gff 1/6 spec!
% 04Nov2009 WCS Modified to integrate with the MATLAB SAR toolbox
%                  Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%======================================================================
% Get a valid file, assure it's at beginning of file, and open with
% the correct "endian" interpretation:
%======================================================================
% Only use the error dialogs if this was the top-level called function!
UseWarnMsgs = length(dbstack) < 2;

if (nargin >= 2)
   
   if isunix,
      if (PathName(length(PathName)) ~= '/') & (PathName(length(PathName)) ~= '\'),
         PathName = [PathName '/'];
      end
      FilterSpec = [[PathName '*.GFF|'],[PathName '*.gff']];
   else
      if (PathName(length(PathName)) ~= '/') & (PathName(length(PathName)) ~= '\'),
         PathName = [PathName '\'];
      end
      FilterSpec = [PathName '*.gff;*.GFF'];
   end
   
   [gffName, gffPath] = uigetfile(FilterSpec, 'Select GFF file to open');
   if (gffName == 0)
      fid_out = -1;
      Header = [];
      gffName = [];
      gffPath = [];
      return
   else
      fid = fopen([gffPath gffName],'r','ieee-be');
   end
elseif (nargin >=1)
   if exist( FileName, 'file' ),
      fid = fopen( FileName );
      if fid < 1,
         
         if UseWarnMsgs,
            errordlg(['GSAT File ' gffPath gffName ' cannot be opened.'],'read_gff_header');
         end
         
         fid_out = fid;
         Header = [];
         gffName = [];
         gffPath = [];
         return
      end
   else
      
      if UseWarnMsgs,
         errordlg(['GSAT File ' FileName ' does not exist.'],'read_gff_header');
      end
      
      fid_out = -2;
      Header = [];
      gffName = [];
      gffPath = [];
      return
   end
   gffName = FileName;
   gffPath = [];
else % nargin = 0
   if isunix,
      FilterSpec = ['*.GFF|','*.gff'];
   else
      FilterSpec = '*.gff;*.GFF';
   end
   
   [gffName, gffPath] = uigetfile(FilterSpec, 'Select GFF file to open');
   if (gffName == 0)
      fid_out = -1;
      Header = [];
      gffName = [];
      gffPath = [];
      return
   else
      fid = fopen([gffPath gffName],'r','ieee-be');
   end
end

% Check that it's a GSAT file:
StructID = char(fread(fid,[1,7],'char'));
if ~strcmp(StructID,'GSATIMG'),
   
   if UseWarnMsgs,
      errordlg('Selected file is not in GFF.','read_gff_header');
   end
   
   fclose(fid);
   fid_out = -3;
   Header = [];
   gffName = [];
   gffPath = [];
   return;
end

% Check for "endian"
fseek(fid,54,'bof');
Header.Endian = fread(fid,1,'ushort');  % 0 = little-endian

% Reopen the file with the proper endian interpretation:
if ( ~Header.Endian )
   fclose(fid);
   fid = fopen([gffPath gffName],'r','ieee-le');
else
   fclose(fid);
   fid = fopen([gffPath gffName],'r','ieee-be');
end


%======================================================================
% The file is now open correctly.  Read the parameters:
%======================================================================

% Check for "endian"
fseek(fid,8,'bof');
Header.Version_Minor     = fread(fid,1,'ushort');
Header.Version_Major     = fread(fid,1,'ushort');
Header.Length            = fread(fid,1,'ulong');
Header.Creator_Length    = fread(fid,1,'ushort');
Header.Creator           = char(fread(fid,[1,Header.Creator_Length],'char'));
                           fread(fid,24-Header.Creator_Length,'char');

if (ftell(fid) ~= 42)
   warndlg('file marker not at 42 to start date & time.');
end

Header.DateTime          = fread(fid,6,'ushort'); % yr, mo, da, hr, min, sec.

Header.Endian            = fread(fid,1,'ushort');  % 0 = little-endian
if (( Header.Version_Major == 1 ) & ( Header.Version_Minor > 7 ) | ...
      (Header.Version_Major > 1)),
       
      Header.BytesPerPixel     = fread(fid,1,'float32'); 
else
      Header.BytesPerPixel     = fread(fid,1,'uint32'); 
end



Header.FrameCnt             = fread(fid,1,'ulong');
Header.ImageType            = fread(fid,1,'ulong');
Header.RowMajor             = fread(fid,1,'ulong');
Header.RgCnt                  = fread(fid,1,'ulong');
Header.AzCnt                  = fread(fid,1,'ulong');
Header.ScaleExponent     = fread(fid,1,'long');
Header.ScaleMantissa     = fread(fid,1,'long');
Header.OffsetExponent    = fread(fid,1,'long');
Header.OffsetMantissa    = fread(fid,1,'long');

                           fread(fid,32,'uchar'); % Throw away "Res2" required filler

if (ftell(fid) ~= 128)
   warndlg('file marker not at 128 to start comment.');
end

Header.Comment_Length    = fread(fid,1,'ushort');
Header.Comment           = char(fread(fid,[1,Header.Comment_Length],'char'));
                           fread(fid,166-Header.Comment_Length,'char');

Header.ImagePlane        = fread(fid,1,'ulong');
Header.RgPixelSz         = fread(fid,1,'ulong') * 2^(-16);
Header.AzPixelSz         = fread(fid,1,'ulong') * 2^(-16);
Header.AzOverlap         = fread(fid,1,'long') * 2^(-16);

Header.SRPLat            = fread(fid,1,'long') * 2^(-23);
Header.SRPLong           = fread(fid,1,'long') * 2^(-23);
Header.SRPAlt            = fread(fid,1,'long') * 2^(-16);

Header.RFOA             = fread(fid,1,'long') * 2^(-23);

Header.XtoSRP            = fread(fid,1,'long') * 2^(-16);

                           fread(fid,32,'uchar'); % Throw away "Res3" required filler

if (ftell(fid) ~= 364)
   warndlg('file marker not at 364 to start phase history ID.');
end

Header.PhaseName_Length  = fread(fid,1,'ushort');
Header.PhaseName         = char(fread(fid,[1,Header.PhaseName_Length],'char'));
                           fread(fid,128-Header.PhaseName_Length,'char');
Header.ImageName_Length  = fread(fid,1,'ushort');
Header.ImageName         = char(fread(fid,[1,Header.ImageName_Length],'char'));
                           fread(fid,128-Header.ImageName_Length,'char');
                           
Header.LookCnt           = fread(fid,1,'ulong');
Header.ParamRefAp        = fread(fid,1,'ulong');
Header.ParamRefPos       = fread(fid,1,'ulong');

Header.GrzAngle          = fread(fid,1,'ulong') * 2^(-23);
Header.Squint            = fread(fid,1,'long') * 2^(-23);
Header.GTA               = fread(fid,1,'long') * 2^(-23);
Header.RgBeamCtr         = fread(fid,1,'ulong') * 2^(-8);
Header.FlightTime        = fread(fid,1,'ulong') * 10^(-3);

if (( Header.Version_Major == 1 ) & ( Header.Version_Minor > 5 ) | ...
      (Header.Version_Major > 1)),
   Header.RgChirpRate       = fread(fid,1,'float');
else
   Header.RgChirpRate       = fread(fid,1,'long') * 2^(+12);
end

Header.XtoStart          = fread(fid,1,'long') * 2^(-16);

Header.MoCompMode        = fread(fid,1,'ulong');

Header.V_X               = fread(fid,1,'ulong') * 2^(-16);

Header.APCLat            = fread(fid,1,'long') * 2^(-23);
Header.APCLong           = fread(fid,1,'long') * 2^(-23);
Header.APCAlt            = fread(fid,1,'ulong') * 2^(-16);

Header.Calparm           = fread(fid,1,'ulong') * 2^(-24);

Header.LogicalBlkAdrr    = fread(fid,1,'ulong');

if (ftell(fid) ~= 692)
   warndlg('file marker not at 692 to start Azimuth Resolution.');
end

Header.AzResolution      = fread(fid,1,'ulong') * 2^(-16);
Header.RgResolution      = fread(fid,1,'ulong') * 2^(-16);

% Even though draft 1-6 says this is unsigned, it is actually SIGNED.
Header.DesSigmaN         = fread(fid,1,'long') * 2^(-23);

Header.DesGrazAngle      = fread(fid,1,'ulong') * 2^(-23);
Header.DesSquint         = fread(fid,1,'long') * 2^(-23);
Header.DesRng            = fread(fid,1,'ulong') * 2^(-8);

Header.SceneTrackAngle   = fread(fid,1,'long') * 2^(-23);

if (ftell(fid) ~= 720)
   warndlg('file marker not at 720 to start User Specified Data.');
end

Header.UserParam         = fread(fid,48,'uchar');  % Skip User Specified Data
                           
if (ftell(fid) ~= 768)
   warndlg('file marker not at 768 to start CCD Parameters.');
end

Header.CoarseSNR         = fread(fid,1,'long');  % CCD Parameters
Header.CoarseAzSubSamp   = fread(fid,1,'long');   
Header.CoarseRngSubSamp  = fread(fid,1,'long');  
Header.MaxAzShift        = fread(fid,1,'long');   
Header.MaxRngShift       = fread(fid,1,'long');   
Header.CoarseDltAz       = fread(fid,1,'long');   
Header.CoarseDltRng      = fread(fid,1,'long');   
Header.TotProcs          = fread(fid,1,'long');   
Header.TptBoxCMode       = fread(fid,1,'long');   
Header.SNRThresh         = fread(fid,1,'long');   
Header.RngSize           = fread(fid,1,'long');  
Header.MapBoxSize        = fread(fid,1,'long');   
Header.BoxSize           = fread(fid,1,'long');   
Header.BoxSpc            = fread(fid,1,'long');  
Header.TotTPts           = fread(fid,1,'long');  
Header.GoodTPts          = fread(fid,1,'long');  
Header.RndSeed           = fread(fid,1,'long');   
Header.RngShift          = fread(fid,1,'long');   
Header.AzShift           = fread(fid,1,'long');   
Header.SumXRamp          = fread(fid,1,'long');  
Header.SumYRamp          = fread(fid,1,'long');  


if (ftell(fid) ~= 852)
   warndlg('file marker not at 852 to start General Fields for RTV.');
end

Header.Cy9kTapeBlock     = fread(fid,1,'ulong');

Header.NominalCenterFreq = fread(fid,1,'float');
                           
Header.ImageFlags        = fread(fid,1,'ulong');

Header.LineNumber        = fread(fid,1,'ulong');
Header.PatchNumber       = fread(fid,1,'ulong');

Header.Lambda0           = fread(fid,1,'float');
Header.SRngPixSpace      = fread(fid,1,'float');
Header.DoppPixSpace      = fread(fid,1,'float');
Header.DoppOffset        = fread(fid,1,'float');
Header.DoppRngScale      = fread(fid,1,'float');
                           
Header.MuxTimeDelay      = fread(fid,1,'float');
                           
Header.APCXECEF          = fread(fid,1,'double');
Header.APCYECEF          = fread(fid,1,'double');
Header.APCZECEF          = fread(fid,1,'double');

Header.VxECEF            = fread(fid,1,'float');
Header.VyECEF            = fread(fid,1,'float');
Header.VzECEF            = fread(fid,1,'float');

if (ftell(fid) ~= 932)
   warndlg('file marker not at 932 to start Phase Calibration for IFSAR.');
end

Header.PhaseCal          = fread(fid,1,'float');

Header.SRPxECEF          = fread(fid,1,'double');
Header.SRPyECEF          = fread(fid,1,'double');
Header.SRPzECEF          = fread(fid,1,'double');

Header.Res5              = fread(fid,64,'uchar');  % filler Res5


%START REV 1.8 CHANGES


if (( Header.Version_Major == 1 ) & ( Header.Version_Minor > 7 ) | ...
      (Header.Version_Major > 1)),

  
Header.HeaderLen1        = fread(fid,1,'ulong');

Header.ImgDate.Year        = fread(fid,1,'ushort');
Header.ImgDate.Month       = fread(fid,1,'ushort');
Header.ImgDate.Day         = fread(fid,1,'ushort');
Header.ImgDate.Hour        = fread(fid,1,'ushort');
Header.ImgDate.Minute      = fread(fid,1,'ushort');
Header.ImgDate.Second      = fread(fid,1,'ushort');

Header.CompFileName      = fread(fid,128,'uchar');
Header.RefFileName       = fread(fid,128,'uchar');
Header.IEPlatform        = fread(fid,24,'uchar');
Header.IEProcID          = fread(fid,12,'uchar');
Header.IERadarModel      = fread(fid,12,'uchar');
Header.IERadarID         = fread(fid,1,'ulong');
Header.IESWID            = fread(fid,24,'uchar');
Header.IFPlatform        = fread(fid,24,'uchar');
Header.IFProcID          = fread(fid,12,'uchar');
Header.IFRadarModel      = fread(fid,12,'uchar');
Header.IFRadarID         = fread(fid,1,'ulong');
Header.IFSWID            = fread(fid,24,'uchar');
Header.IFAlgo            = fread(fid,8,'uchar');
Header.PHPlatform        = fread(fid,24,'uchar');
Header.PHProcID          = fread(fid,12,'uchar');
Header.PHRadarModel      = fread(fid,12,'uchar');
Header.PHRadarID         = fread(fid,1,'ulong');
Header.PHSWID            = fread(fid,24,'uchar');

Header.PHDataRcd         = fread(fid,1,'ulong');

Header.ProcProduct       = fread(fid,1,'ulong');
Header.MissionText       = fread(fid,8,'uchar');

Header.PHSource          = fread(fid,1,'ulong');

Header.GPSWeek           = fread(fid,1,'ulong');
Header.DataCollectReqH   = fread(fid,14,'uchar');
Header.Res6              = fread(fid,2,'uchar');
Header.GridName          = fread(fid,24,'uchar');

Header.PixValLinearity   = fread(fid,1,'ulong');

Header.ComplexOrReal     = fread(fid,1,'ulong');

Header.BitsPerMagnitude  = fread(fid,1,'ushort');
Header.BitsPerPhase      = fread(fid,1,'ushort');

Header.ComplexOrderType  = fread(fid,1,'ulong');

Header.PixDataType       = fread(fid,1,'ulong');

Header.ImageLength       = fread(fid,1,'ulong');

Header.ImageCmpScheme    = fread(fid,1,'ulong');

Header.APBO              = fread(fid,1,'float');
Header.AsaPitch          = fread(fid,1,'float');
Header.AsaSquint         = fread(fid,1,'float');
Header.DsaPitch          = fread(fid,1,'float');
Header.IRA               = fread(fid,1,'float');
Header.RxPolarization    = fread(fid,2,'float');
Header.TxPolarization    = fread(fid,2,'float');
Header.VxAvg             = fread(fid,1,'float');
Header.VyAvg             = fread(fid,1,'float');
Header.VzAvg             = fread(fid,1,'float');
Header.APCxAvg           = fread(fid,1,'float');
Header.APCyAvg           = fread(fid,1,'float');
Header.APCzAvg           = fread(fid,1,'float');
Header.AveragingTime     = fread(fid,1,'float');
Header.Dgta              = fread(fid,1,'float');
Header.VelocY            = fread(fid,1,'ulong')*2^(-16);
Header.VelocZ            = fread(fid,1,'ulong')*2^(-16);
Header.Ba                = fread(fid,1,'float');
Header.Be                = fread(fid,1,'float');

Header.AzGeomCorr        = fread(fid,1,'ulong');

Header.RngGeomCorr       = fread(fid,1,'ulong');

Header.AzWinFacBW        = fread(fid,1,'float');
Header.RngWinFacBW       = fread(fid,1,'float');
Header.AzWinID           = fread(fid,48,'uchar');
Header.RngWinID          = fread(fid,48,'uchar');
Header.KeepOutViolPrcnt  = fread(fid,1,'float');
Header.AzCoeff           = fread(fid,6,'float');
Header.PosUncertDown     = fread(fid,1,'float');
Header.PosUncertE        = fread(fid,1,'float');
Header.PosUncertN        = fread(fid,1,'float');

Header.NavAidingType     = fread(fid,1,'ulong');

Header.TwoDNLPhaseCoeffs = fread(fid,10,'float');
Header.ClutterSNRThresh  = fread(fid,1,'float');
Header.ElevationCoeff    = fread(fid,9,'float');
Header.MonopulseCoeff    = fread(fid,12,'float');
Header.TwistPntErrPrcnt  = fread(fid,1,'float');
Header.TiltPntErrPrcnt   = fread(fid,1,'float');
Header.AzPntErrPrcnt     = fread(fid,1,'float');
Header.SigmaN            = fread(fid,1,'ulong')*2^(-23);
Header.TakeNum           = fread(fid,1,'ulong');
Header.IFSARFlags        = fread(fid,5,'ulong');
Header.MuThreshold       = fread(fid,1,'float');

Header.GffAppType        = fread(fid,1,'ulong');

Header.Res7              = fread(fid,8,'uchar');
end


 if (ftell(fid) ~= Header.Length)
    warndlg('File marker not at header length - moving to end.');
   % Position file at end of header:
    fseek(fid,Header.Length,'bof');
end

if (nargout > 3)
   fid_out = fid;
else
   fclose(fid);
end

return

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////