function [ meta ] = read_nitf_meta( filename )
%READ_NITF_META Read metadata from NITF into a structure
%
%   NITFHDR = READ_NITF_META(NITFFILE) calls functions to parse NITF file
%   segments.
%
% Written by: Matt Donath, NGA, matthew.b.donath@nga.ic.gov
% Modified by: Wade Schwartzkopf, NGA/IDT
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

%% OPEN FILE
fid=fopen(filename,'r','b'); % NITF always big-endian

%% READ NITF FILE HEADER
meta.filehdr = read_nitf_filehdr(fid);

%% READ IMAGE SEGMENTS
nextImgSubhdroffset=meta.filehdr.HL;
for i = 1:meta.filehdr.NUMI
    fseek(fid,nextImgSubhdroffset,'bof'); % Jump to ith image segment
    meta.imagesubhdr{i}=read_nitf_imgsubhdr(fid);
    nextImgSubhdroffset=nextImgSubhdroffset+meta.filehdr.LISH(i)+meta.filehdr.LI(i);
end
fseek(fid,nextImgSubhdroffset,'bof'); % Jump to end of image segments

%% SKIP GRAPHIC SEGMENTS
if isfield(meta.filehdr,'NUMS') && meta.filehdr.NUMS > 0
    for seg = 1:meta.filehdr.NUMS
        fseek(fid, meta.filehdr.LSSH(seg), 'cof'); % Skip graphic segment subheaders
        fseek(fid, meta.filehdr.LS(seg), 'cof'); % Skip graphic segment data
    end
end

%% SKIP RESERVED SEGMENTS
if isfield(meta.filehdr,'NUMX') && meta.filehdr.NUMX > 0
    for seg = 1:meta.filehdr.NUMX
        fseek(fid, meta.filehdr.LXSH(seg), 'cof'); % Skip reserved segment subheaders
        fseek(fid, meta.filehdr.LX(seg), 'cof'); % Skip reserved segment data
    end
end

%% SKIP LABEL SEGMENTS
if isfield(meta.filehdr,'NUML') && meta.filehdr.NUML > 0
    for seg = 1:meta.filehdr.NUML
        fseek(fid, meta.filehdr.LLSH(seg), 'cof'); % Skip label segment subheaders
        fseek(fid, meta.filehdr.LL(seg), 'cof'); % Skip label segment data
    end
end

%% SKIP TEXT SEGMENTS
if isfield(meta.filehdr,'NUMT') && meta.filehdr.NUMT > 0
    for seg = 1:meta.filehdr.NUMT
%         % Skip text segments
%         fseek(fid, meta.filehdr.LTSH(seg), 'cof'); % Skip text segment subheaders
%         fseek(fid, meta.filehdr.LT(seg), 'cof'); % Skip text segment data
        % Read text segments
        meta.txtsubhdr{seg}=read_nitf_txtsubhdr(fid);
        meta.txtseg{seg} = fread(fid,meta.filehdr.LT(seg),'uint8=>char')';
    end
end

%% READ DATA EXTENSION SEGMENT SUBHEADERS
curr_pos = ftell(fid);
if isfield(meta.filehdr,'NUMDES') && meta.filehdr.NUMDES > 0
    for seg = 1:meta.filehdr.NUMDES
        fseek(fid, curr_pos, 'bof');
        fld = strcat('DESSubhdr',num2str(seg));
        meta = readDESegSubhdr(meta,fid,fld);
        curr_pos = curr_pos + meta.filehdr.LDSH(seg) + meta.filehdr.LD(seg);
    end
end

%% READ RESERVED EXTENSION SUBHEADERS
if isfield(meta.filehdr,'NUMRES') && meta.filehdr.NUMRES > 0
    for seg = 1:meta.filehdr.NUMRES
        fld = strcat('RESSubhdr',num2str(seg));
        meta = readRESegSubhdr(meta,fid,fld);
    end
end

%% Close file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////