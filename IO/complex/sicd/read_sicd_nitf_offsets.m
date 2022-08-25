function [ NITF_meta ] = read_sicd_nitf_offsets( filename )
%READ_SICD_NITF_OFFSETS Read metadata from Sensor Independent Complex Data
%(SICD) file, version 0.3 and above
%
% SICD (version >= 0.3) is stored in a NITF container.  NITF is a
% complicated format that involves lots of fields and configurations
% possibilities. Fortunately, SICD only really uses a small, specific
% portion of the NITF format.  This function extracts only the few parts of
% the NITF metadata necessary for reading a SICD NITF file.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Open file
fid = fopen(filename);

%% Read NITF file header
if ~strcmp(fread(fid,9,'uint8=>char')','NITF02.10') % Check format
    fclose(fid);
    error('READ_SICD_META:INVALID_FILE_FORMAT','Only NITF version 2.1 recognized.');
end
fseek(fid,354,'bof'); % Offset to first field of interest
HL = str2double(fread(fid,6,'uint8=>char')'); % File header length
NUMI = str2double(fread(fid,3,'uint8=>char')'); % Number of image segments
if NUMI > 0
    ImgSubhdrLngths=zeros(NUMI,1); % Length of the image segment subheader
    ImgLngths=zeros(NUMI,1); % Length of the image segment data
    NITF_meta.minimal.imageSegmentOffsets=zeros(NUMI,1); % Offset to image data from beginning of file (in bytes)
    NITF_meta.minimal.imageSegmentRows=zeros(NUMI,1); % Number of rows in each image segment (in case data is spread across multiple image segments)
    NITF_meta.minimal.imageSegmentColumns=zeros(NUMI,1); % Number of rows in each image segment (in case data is spread across multiple image segments)
    for i = 1:NUMI
        ImgSubhdrLngths(i) = str2double(fread(fid,6,'uint8=>char')');
        NITF_meta.minimal.imageSegmentOffsets(i)=HL+sum(ImgSubhdrLngths)+sum(ImgLngths);
        ImgLngths(i) = str2double(fread(fid,10,'uint8=>char')');
    end
else
    NITF_meta.minimal.imageSegmentOffsets=[];
    NITF_meta.minimal.imageSegmentRows=[];
    NITF_meta.minimal.imageSegmentColumns=[];
end
for i=1:2 % Jump over segments not used by SICD to DES description
    segment_length = str2double(fread(fid,3,'uint8=>char')');
    if(isfinite(segment_length)&&(segment_length>0))
        error('READ_SICD_META:INVALID_NITF_SEGMENT','SICD does not allow for graphics or reserved extension segments.');
    end
end
% Jump over text segments
NUMT = str2double(fread(fid,3,'uint8=>char')');
LTSH = zeros(NUMT,1); LT = zeros(NUMT,1);
for lp = 1:NUMT
    LTSH(lp) = str2double(fread(fid,4,'uint8=>char')'); % Length of text segment subheader
    LT(lp) = str2double(fread(fid,5,'uint8=>char')'); % Length of text segment data
end
NUMDES = str2double(fread(fid,3,'uint8=>char')');  % Number of data extension segments
LDSH = zeros(NUMDES,1); LD = zeros(NUMDES,1);
for lp = 1:NUMDES
    LDSH(lp) = str2double(fread(fid,4,'uint8=>char')'); % Length of data extension segment subheader
    LD(lp) = str2double(fread(fid,9,'uint8=>char')'); % Length of data extension segment data
end
NITF_meta.minimal.desLengths = LD;
NITF_meta.minimal.desOffsets = HL + sum(ImgSubhdrLngths) + sum(ImgLngths) + ...
    sum(LTSH) + sum(LT) + cumsum(LDSH) + [0; cumsum(LD(1:(end-1)))];

%% Get number of rows for each image segment from image segment headers
nextImgSubhdroffset=HL;
for i = 1:NUMI
    fseek(fid,nextImgSubhdroffset,'bof'); % Jump to ith image segment
    fseek(fid,333,'cof'); % Jump to number of rows field
    NITF_meta.minimal.imageSegmentRows(i)=str2double(fread(fid,8,'uint8=>char')');
    NITF_meta.minimal.imageSegmentColumns(i)=str2double(fread(fid,8,'uint8=>char')');
    nextImgSubhdroffset=nextImgSubhdroffset+ImgSubhdrLngths(i)+ImgLngths(i);
end

%% Close file
fclose(fid);

% If you want the full NITF header data, you can use the NITF functions
% from any other NITF toolbox, but its not really necessary for anything in
% SICD.
% NITF_meta.IPT = nitfinfo(filename); % MATLAB Image Processing Toolbox
% NITF_meta.matlab_sar = read_nitf_meta(filename); % NITF reading in the MATLAB SAR Toolbox

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////