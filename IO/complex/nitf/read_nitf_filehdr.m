function [filehdr] = read_nitf_filehdr(fid)
%READ_NITF_FILEHDR Read and populate fields into file header structure.
%   [FILEHDR] = READ_NITF_FILEHDR(FID) read the nitf file header and 
%   place the metadata into a structure. Handles NITF versions 2.0 and 2.1.
%   Takes a file identifier pointing to the beginning of the file header
%   (which is the beginning of the file). 
%
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

%% COMMON FIELDS

filehdr.FHDR = fread(fid,9,'uint8=>char')';
filehdr.CLEVEL = fread(fid,2,'uint8=>char')';
filehdr.STYPE = strtrim(fread(fid,4,'uint8=>char')');
filehdr.OSTAID = strtrim(fread(fid,10,'uint8=>char')');
filehdr.FDT = fread(fid,14,'uint8=>char')';
filehdr.FTITLE = strtrim(fread(fid,80,'uint8=>char')');
filehdr.FSCLAS = fread(fid,1,'uint8=>char')';

%% VERSION 2.1

if strcmp(filehdr.FHDR,'NITF02.10')
    
    % Security info
    filehdr.FSCLSY = strtrim(fread(fid,2,'uint8=>char')');
    filehdr.FSCODE = strtrim(fread(fid,11,'uint8=>char')');
    filehdr.FSCTLH = strtrim(fread(fid,2,'uint8=>char')');
    filehdr.FSREL = strtrim(fread(fid,20,'uint8=>char')');
    filehdr.FSDCTP = strtrim(fread(fid,2,'uint8=>char')');
    filehdr.FSDCDT = strtrim(fread(fid,8,'uint8=>char')');
    filehdr.FSDCXM = strtrim(fread(fid,4,'uint8=>char')');
    filehdr.FSDG = strtrim(fread(fid,1,'uint8=>char')');
    filehdr.FSDGT = strtrim(fread(fid,8,'uint8=>char')');
    filehdr.FSCLTX = strtrim(fread(fid,43,'uint8=>char')');
    filehdr.FSCATP = strtrim(fread(fid,1,'uint8=>char')');
    filehdr.FSCAUT = strtrim(fread(fid,40,'uint8=>char')');
    filehdr.FSCRSN = strtrim(fread(fid,1,'uint8=>char')');
    filehdr.FSSRDT = strtrim(fread(fid,8,'uint8=>char')');
    filehdr.FSCTLN = strtrim(fread(fid,15,'uint8=>char')');
    filehdr.FSCOP = fread(fid,5,'uint8=>char')';
    filehdr.FSCPYS = fread(fid,5,'uint8=>char')';
    filehdr.ENCRYP = fread(fid,1,'uint8=>char')';
    filehdr.FBKGC = uint8(fread(fid,3,'uint8=>char')');
    filehdr.ONAME = strtrim(fread(fid,24,'uint8=>char')');
    filehdr.OPHONE = strtrim(fread(fid,18,'uint8=>char')');
    
    % File and header lengths
    filehdr.FL = str2double(fread(fid,12,'uint8=>char')'); 
    filehdr.HL = str2double(fread(fid,6,'uint8=>char')');
    
    % Image Segments
    filehdr.NUMI = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMI > 0
        for lp = 1:filehdr.NUMI
            filehdr.LISH(lp) = str2double(fread(fid,6,'uint8=>char')');
            filehdr.LI(lp) = str2double(fread(fid,10,'uint8=>char')');
        end
    end
    
    % Graphic segments
    filehdr.NUMS = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMS > 0
        for lp = 1:filehdr.NUMS
            filehdr.LSSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LS(lp) = str2double(fread(fid,6,'uint8=>char')');
        end
    end
    
    % Reserved extension segments
    filehdr.NUMX = str2double(fread(fid,3,'uint8=>char')');
    
    % Text Segments
    filehdr.NUMT = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMT > 0
        for lp = 1:filehdr.NUMT
            filehdr.LTSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LT(lp) = str2double(fread(fid,5,'uint8=>char')');
        end
    end

    % Data Extension segments
    filehdr.NUMDES = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMDES > 0
        for lp = 1:filehdr.NUMDES
            filehdr.LDSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LD(lp) = str2double(fread(fid,9,'uint8=>char')');
        end
    end

    % Reserved Extension segments
    filehdr.NUMRES = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMRES > 0
        for lp = 1:filehdr.NUMRES
            filehdr.LRESH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LRE(lp) = str2double(fread(fid,7,'uint8=>char')');
        end
    end
    
    % User defined headers
    filehdr.UDHDL = str2double(fread(fid,5,'uint8=>char')');
    if filehdr.UDHDL > 0
        filehdr.UDHOFL = str2double(fread(fid,3,'uint8=>char')');
        filehdr.UDHD = str2double(fread(fid,filehdr.UDHDL-3,'uint8=>char')');
    else
        filehdr.UDHOFL = 0;
        filehdr.UDHD = 0;
    end
    
    % Extended headers
    filehdr.XHDL = str2double(fread(fid,5,'uint8=>char')');
    if filehdr.XHDL > 0
        filehdr.XHDLOFL = str2double(fread(fid,3,'uint8=>char'))';
        filehdr.XHD = read_ext_hdr(fid,filehdr.XHDL);
    else
        filehdr.XHDLOFL = 0;
        filehdr.XHD = 0;
    end

%% VERSION 2.0

elseif strcmp(filehdr.FHDR,'NITF02.00')
    % Security info
    filehdr.FSCODE = strtrim(fread(fid,40,'uint8=>char')');
    filehdr.FSCTLH = strtrim(fread(fid,40,'uint8=>char')');
    filehdr.FSREL =  strtrim(fread(fid,40,'uint8=>char')');
    filehdr.FSCAUT = strtrim(fread(fid,20,'uint8=>char')');
    filehdr.FSCTLN = strtrim(fread(fid,20,'uint8=>char')');
    filehdr.FSDWNG = strtrim(fread(fid,6,'uint8=>char')');
    
    if strcmp(filehdr.FSDWNG,'999998')
        filehdr.FSDEVT = strtrim(fread(fid,40,'uint8=>char')');
    
    end
    filehdr.FSCOP =  strtrim(fread(fid,5,'uint8=>char')');
    filehdr.FSCPYS = strtrim(fread(fid,5,'uint8=>char')');
    filehdr.ENCRYP = strtrim(fread(fid,1,'uint8=>char')');
    filehdr.FBKGC = fread(fid,3,'uint8')';
    filehdr.ONAME =  strtrim(fread(fid,24,'uint8=>char')');
    filehdr.OPHONE = strtrim(fread(fid,18,'uint8=>char')');
    
    % File and header lengths
    filehdr.FL = str2double(fread(fid,12,'uint8=>char')');
    filehdr.HL = str2double(fread(fid,6,'uint8=>char')');

    % Image Segments
    filehdr.NUMI = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMI > 0
        for lp = 1:filehdr.NUMI
            filehdr.LISH(lp) = str2double(fread(fid,6,'uint8=>char')');
            filehdr.LI(lp) = str2double(fread(fid,10,'uint8=>char')');
        end
    end

    % Symbol segments
    filehdr.NUMS = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMS > 0
        for lp = 1:filehdr.NUMS
            filehdr.LSSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LS(lp) = str2double(fread(fid,6,'uint8=>char')');
        end
    end
    
    % Label segments
    filehdr.NUML = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUML > 0
        for lp = 1:filehdr.NUML
            filehdr.LLSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LL(lp) = str2double(fread(fid,3,'uint8=>char')');
        end
    end

    % Text segments
    filehdr.NUMT = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMT > 0
        for lp = 1:filehdr.NUMT
            filehdr.LTSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LT(lp) = str2double(fread(fid,5,'uint8=>char')');
        end
    end

    % Data Extension segments
    filehdr.NUMDES = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMDES > 0
        for lp = 1:filehdr.NUMDES
            filehdr.LDSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LD(lp) = str2double(fread(fid,9,'uint8=>char')');
        end
    end
    
    % Reserved Extension segments
    filehdr.NUMRES = str2double(fread(fid,3,'uint8=>char')');
    if filehdr.NUMRES > 0
        for lp = 1:filehdr.NUMRES
            filehdr.LRSH(lp) = str2double(fread(fid,4,'uint8=>char')');
            filehdr.LR(lp) = str2double(fread(fid,7,'uint8=>char')');
        end
    end

    % User defined headers
    filehdr.UDHDL = str2double(fread(fid,5,'uint8=>char')');
    if filehdr.UDHDL > 0
        filehdr.UDHOFL = fread(fid,3,'uint8=>char')';
        filehdr.UDHD = fread(fid,filehdr.UDHDL-3,'uint8=>char')';
        
    end
    
    % Extended headers
    filehdr.XHDL = str2double(fread(fid,5,'uint8=>char')');
    if filehdr.XHDL > 0
       	filehdr.XDHOFL = fread(fid,3,'uint8=>char')';
        filehdr.XHD = fread(fid,filehdr.XHDL-3,'uint8=>char')';
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////