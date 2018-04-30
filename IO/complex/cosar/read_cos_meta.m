function [ meta ] = read_cos_meta( filename )
%READ_COS_META Read metadata from TerraSAR-X COSAR file format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open file
fid=fopen(filename,'r','b');

% Read metadata
av_format='int32'; % Annotation value binary format
meta.bib=fread(fid,1,av_format); % The number of bytes in the actual burst
%   (BIB = Bytes In Burst). Including the annotation and valid only for
%   ScanSAR bursts.
meta.rsri=fread(fid,1,av_format); % A range index, giving the relative range
%   location on a virtual common raster with the ADC sampling (its rate is
%   approx. 330MHz) of the bursts first range sample with respect to the
%   reference value (RSRI = Range Sample Relative Index).
meta.rs=fread(fid,1,av_format); % The length of a range line given in
%   samples.  This value has to be same for all bursts and is repeated at
%   every burst (RS = Range Samples).
meta.az=fread(fid,1,av_format); % The length of an azimuth column of the
%   actual burst given in samples. This value may vary from burst to burst
%   (AS = Azimuth Samples).
meta.bi=fread(fid,1,av_format); % The index number of the burst (BI = Burst
%   Index).
meta.rtnb=fread(fid,1,av_format); % The total number of bytes in a line in
%   range direction (the “width” of the entire file including the annotation
%   bytes). As the TNL, this parameter is given only once in the first line
%   of the file (RTNB = Rangeline Total Number of Bytes).
meta.tnl=fread(fid,1,av_format); % The extent in azimuth direction (the
%   “height” of the entire file including the annotation lines).  This
%   parameter is given only once in the first line of the file in order to
%   facilitate the reading of the file and replaced by the special filler
%   value for the other bursts (TNL = Total Number of Lines). The file size
%   can be derived from RTNB times TNL.
meta.format=fread(fid,4,'uchar=>char');
meta.version=fread(fid,1,av_format); % For the convenience of
%   multi-format reader software the following 2 samples identify the file
%   format (not visible in Figure 4-4). The first sample reads hex.
%   43534152 which is the ASCII string CSAR and the second sample gives a
%   version number.
meta.oversample=fread(fid,1,av_format); % The following sample gives the
%   oversampling factor of the RSRI sample position with respect to the
%   current range sampling (1 for 330MHz, 2 for 165MHz or 3 for 110MHz).
meta.scalingrate=fread(fid,1,'float64'); % The following two samples
%   contain the 8-byte floating point value (MSB order) of the inverse 
%   SPECAN scaling rate 1/k applied in processing of the burst. This
%   information may facilitate interferometric processing but it is not
%   meaningfull for Stripmap modes (1/k -> 0).

% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////