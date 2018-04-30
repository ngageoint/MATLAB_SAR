function nitfHdr = readRESegSubhdr(nitfHdr,fid,numREs)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
%READIMGSUBHDR Create structure for reserved extension metadata.
%
% >>>> DUE TO LACK OF APPROPRIATE DATA, THIS FUNCTION IS UNTESTED. <<<<

% Written by Matt Donath, NGA, matthew.b.donath@nga.ic.gov
% References:
%  > MIL-STD-2500A, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.0
%  > MIL-STD-2500C, NATIONAL IMAGERY TRANSMISSION FORMAT VERSION 2.1
%  > STDI-0001,     NATIONAL SUPPORT DATA EXTENSIONS (SDE)(VERSION 1.3/CN1)
%                   FOR THE NATIONAL IMAGERY TRANSMISSION FORMAT (NITF)
%  > STDI-0002,     THE COMPENDIUM OF CONTROLLED EXTENSIONS (CE) FOR THE 
%                   NATIONAL IMAGERY TRANSMISSION FORMAT (NITF) VERSION 2.1
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% DETERMINE NITF VERSION
nitfVer = nitfHdr.filehdr.FHDR(6:8);

%% VERSION 2.1

if strcmp(nitfVer ,'2.1')
    nitfHdr.(numREs).RE = fread(fid,2,'uint8=>char')';
    nitfHdr.(numREs).RESID = strtrim(fread(fid,25,'uint8=>char')');

    nitfHdr.(numREs).RESVER = fread(fid,2,'uint8=>char')';
    nitfHdr.(numREs).RECLAS = fread(fid,1,'uint8=>char')';
    nitfHdr.(numREs).RECLSY = fread(fid,2,'uint8=>char')';
    nitfHdr.(numREs).RECODE = strtrim(fread(fid,11,'uint8=>char')');
    nitfHdr.(numREs).RECTLH = strtrim(fread(fid,2,'uint8=>char')');
    nitfHdr.(numREs).REREL =  strtrim(fread(fid,20,'uint8=>char')');
    nitfHdr.(numREs).REDCTP =  strtrim(fread(fid,2,'uint8=>char')');
    nitfHdr.(numREs).REDCDT =  strtrim(fread(fid,8,'uint8=>char')');
    nitfHdr.(numREs).REDCXM =  strtrim(fread(fid,4,'uint8=>char')');
    nitfHdr.(numREs).REDG =  fread(fid,1,'uint8=>char')';
    nitfHdr.(numREs).REDGDT =  strtrim(fread(fid,8,'uint8=>char')');
    nitfHdr.(numREs).RECLTX =  strtrim(fread(fid,43,'uint8=>char')');
    nitfHdr.(numREs).RECATP =  fread(fid,1,'uint8=>char')';
    nitfHdr.(numREs).RECAUT = strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numREs).RECRSN = fread(fid,1,'uint8=>char')';
    nitfHdr.(numREs).RESRDT = strtrim(fread(fid,8,'uint8=>char')');
    nitfHdr.(numREs).RECTLN = strtrim(fread(fid,15,'uint8=>char')');
        
    nitfHdr.(numREs).RESSHL = str2double(fread(fid,4,'uint8=>char')');
    if nitfHdr.(numREs).RESSHL > 0
        nitfHdr.(numREs).RESSHF = fread(fid,nitfHdr.(numREs).RESSHL,'uint8=>char')';
    end
    
    if strcmp(nitfHdr.(numREs).RESID,'TRE_OVERFLOW')
        TRE = fread(fid,6,'uint8=>char')';
        fseek(fid, -6, 'cof');
        nitfHdr.(numREs).RESDATA = eval(strcat('read',TRE,'(fid)'));
    end

%% VERSION 2.0

elseif strcmp(nitfVer, '2.0')
    nitfHdr.(numREs).RE = fread(fid,2,'uint8=>char')';
    nitfHdr.(numREs).RESTAG = strtrim(fread(fid,25,'uint8=>char')');
    nitfHdr.(numREs).RESVER = fread(fid,2,'uint8=>char')';
    
    % nitfHdr.(numREs).RESSG
    nitfHdr.(numREs).RESSG.RESCLAS = fread(fid,1,'uint8=>char')';
    nitfHdr.(numREs).RESSG.RESCODE = strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numREs).RESSG.RESCTLH = strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numREs).RESSG.RESREL =  strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numREs).RESSG.RESCAUT = strtrim(fread(fid,20,'uint8=>char')');
    nitfHdr.(numREs).RESSG.RESCTLN = strtrim(fread(fid,20,'uint8=>char')');
    nitfHdr.(numREs).RESSG.RESDWNG = strtrim(fread(fid,6,'uint8=>char')');
    if strcmp(nitfHdr.(numREs).RESSG.RESDWNG,'999998')
        nitfHdr.(numREs).RESSG.RESDEVT = strtrim(fread(fid,40,'uint8=>char')');
    end
    
    if strcmp(nitfHdr.(numREs).RESTAG,'Registered Extensions') ||...
       strcmp(nitfHdr.(numREs).RESTAG,'Controlled Extensions')
          nitfHdr.(numREs).RESOFLW = fread(fid,6,'uint8=>char')';
          nitfHdr.(numREs).RESITEM = fread(fid,3,'uint8=>char')';
    end
    
    nitfHdr.(numREs).RESSHL = str2double(fread(fid,4,'uint8=>char')');
    if nitfHdr.(numREs).RESSHL > 0
        nitfHdr.(numREs).RESSHF = fread(fid,nitfHdr.(numREs).RESSHL,'uint8=>char')';
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////