function nitfHdr = readDESegSubhdr(nitfHdr,fid,numDEs)
%
% <<<<<<<<<< CLASSIFICATION: UNCLASSIFIED >>>>>>>>>> 
%
%READDESEGSUBHDR Read and parse data extension segments

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
    nitfHdr.(numDEs).DE = fread(fid,2,'uint8=>char')';
    nitfHdr.(numDEs).DESID = strtrim(fread(fid,25,'uint8=>char')');
    nitfHdr.(numDEs).DESVER = fread(fid,2,'uint8=>char')';
    nitfHdr.(numDEs).DECLAS = fread(fid,1,'uint8=>char')';
    nitfHdr.(numDEs).DESCLSY = fread(fid,2,'uint8=>char')';
    nitfHdr.(numDEs).DESCODE = strtrim(fread(fid,11,'uint8=>char')');
    nitfHdr.(numDEs).DESCTLH = strtrim(fread(fid,2,'uint8=>char')');
    nitfHdr.(numDEs).DESREL =  strtrim(fread(fid,20,'uint8=>char')');
    nitfHdr.(numDEs).DESDCTP =  strtrim(fread(fid,2,'uint8=>char')');
    nitfHdr.(numDEs).DESDCDT =  strtrim(fread(fid,8,'uint8=>char')');
    nitfHdr.(numDEs).DESDCXM =  strtrim(fread(fid,4,'uint8=>char')');
    nitfHdr.(numDEs).DESDG =  fread(fid,1,'uint8=>char')';
    nitfHdr.(numDEs).DESDGDT =  strtrim(fread(fid,8,'uint8=>char')');
    nitfHdr.(numDEs).DESCLTX =  strtrim(fread(fid,43,'uint8=>char')');
    nitfHdr.(numDEs).DESCATP =  fread(fid,1,'uint8=>char')';
    nitfHdr.(numDEs).DESCAUT = strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numDEs).DESCRSN = fread(fid,1,'uint8=>char')';
    nitfHdr.(numDEs).DESSRDT = strtrim(fread(fid,8,'uint8=>char')');
    nitfHdr.(numDEs).DESCTLN = strtrim(fread(fid,15,'uint8=>char')');
    
    if strcmp(nitfHdr.(numDEs).DESID,'TRE_OVERFLOW')
        nitfHdr.(numDEs).DESOFLW = fread(fid,6,'uint8=>char')';
        nitfHdr.(numDEs).DESITEM = str2double(fread(fid,3,'uint8=>char')');
    end
    
    nitfHdr.(numDEs).DESSHL = str2double(fread(fid,4,'uint8=>char')');
    
    if strcmp(nitfHdr.(numDEs).DESID,'TRE_OVERFLOW')
        if nitfHdr.(numDEs).DESSHL > 0
            nitfHdr.(numDEs).DESSHF = fread(fid,nitfHdr.(numDEs).DESSHL,'uint8=>char')';
        end
        TRE = fread(fid,6,'uint8=>char')';
        fseek(fid, -6, 'cof');
        % Check to see if reader function exists for this TRE.
        if exist(strcat('read',TRE,'.m'),'file')
            nitfHdr.(numDEs) = eval(strcat('read',TRE,'(fid)'));
        else
            % Read undefined data
            nitfHdr.(numDEs) = readUndefinedTRE(fid);
        end
    elseif strcmp(nitfHdr.(numDEs).DESID,'STREAMING_FILE_HEADER')
        nitfHdr.(numDEs).SFH_L1 = str2double(fread(fid,7,'uint8=>char')');
        nitfHdr.(numDEs).SFH_DELIM1 = fread(fid,4,'uint8=>char')';
        nitfHdr.(numDEs).SFH_DR = str2double(fread(fid,nitfHdr.(numDEs).SFH_L1,'uint8=>char')');
        nitfHdr.(numDEs).SFH_DELIM2 = fread(fid,4,'uint8=>char')';
        nitfHdr.(numDEs).SFH_L2 = str2double(fread(fid,7,'uint8=>char')');
    elseif strcmp(nitfHdr.(numDEs).DESID,'SICD_XML')
        SICD_xml_string = fread(fid,nitfHdr.filehdr.LD(str2double(numDEs(10:end))),'uint8=>char')';
        nitfHdr.(numDEs).DESDATA = sicdxml2struct( xmlread( java.io.StringBufferInputStream( SICD_xml_string ) ) );
    else % Skip unrecognized DES types
        nitfHdr.(numDEs).DESUSH_OFFSET = ftell(fid); % Record offset, so we can parse user subheaders later
        fseek(fid, nitfHdr.(numDEs).DESSHL, 'cof');
        fseek(fid, nitfHdr.filehdr.LD(str2double(numDEs(10:end))), 'cof');
    end
    
%% VERSION 2.0

elseif strcmp(nitfVer, '2.0')
    nitfHdr.(numDEs).DE = fread(fid,2,'uint8=>char')';
    nitfHdr.(numDEs).DESTAG = strtrim(fread(fid,25,'uint8=>char')');
    nitfHdr.(numDEs).DESVER = fread(fid,2,'uint8=>char')';
    
    % nitfHdr.(numDEs).DESSG
    nitfHdr.(numDEs).DESSG.DESCLAS = fread(fid,1,'uint8=>char')';
    nitfHdr.(numDEs).DESSG.DESCODE = strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numDEs).DESSG.DESCTLH = strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numDEs).DESSG.DESREL =  strtrim(fread(fid,40,'uint8=>char')');
    nitfHdr.(numDEs).DESSG.DESCAUT = strtrim(fread(fid,20,'uint8=>char')');
    nitfHdr.(numDEs).DESSG.DESCTLN = strtrim(fread(fid,20,'uint8=>char')');
    nitfHdr.(numDEs).DESSG.DESDWNG = strtrim(fread(fid,6,'uint8=>char')');
    if strcmp(nitfHdr.(numDEs).DESSG.DESDWNG,'999998')
        nitfHdr.(numDEs).DESSG.DESDEVT = strtrim(fread(fid,40,'uint8=>char')');
    end
    
    if strcmpi(nitfHdr.(numDEs).DESTAG,'Registered Extensions') ||...
            strcmpi(nitfHdr.(numDEs).DESTAG,'Controlled Extensions')
        nitfHdr.(numDEs).DESOFLW = strtrim(fread(fid,6,'uint8=>char')');
        nitfHdr.(numDEs).DESITEM = fread(fid,3,'uint8=>char')';
    end


    nitfHdr.(numDEs).DESSHL = str2double(fread(fid,4,'uint8=>char')');
    if nitfHdr.(numDEs).DESSHL > 0
        nitfHdr.(numDEs).DESSHF = fread(fid,nitfHdr.(numDEs).DESSHL,'uint8=>char')';
    end
    
    if strcmpi(nitfHdr.(numDEs).DESTAG,'Registered Extensions') ||...
            strcmpi(nitfHdr.(numDEs).DESTAG,'Controlled Extensions')
        % READ TRE
        TRE = fread(fid,6,'uint8=>char')';
        fseek(fid, -6, 'cof');

        % Check to see if reader function exists for this TRE.
        if exist(strcat('read',TRE,'.m'),'file')
            nitfHdr.(numDEs) = eval(strcat('read',TRE,'(fid)'));
        else
            % Read undefined data
            nitfHdr.(numDEs) = readUndefinedTRE(fid);
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////