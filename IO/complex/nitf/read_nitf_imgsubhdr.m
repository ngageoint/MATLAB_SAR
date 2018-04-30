function imgsubhdr = read_nitf_imgsubhdr(fid,nitfver)
%READ_NITF_IMGSUBHDR Create structure for required image subheader metadata
%fields.

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

%% DETERMINE NITF VERSION
if nargin<2
    nitfver='2.1'; % Assume version 2.1 unless told otherwise
end

%% VERSION 2.1

if strcmp(nitfver,'2.1')
    imgsubhdr.IM = fread(fid,2,'uint8=>char')';
    imgsubhdr.IID1 = strtrim(fread(fid,10,'uint8=>char')');
    imgsubhdr.IDATIM = fread(fid,14,'uint8=>char')';
    imgsubhdr.TGTID = strtrim(fread(fid,17,'uint8=>char')');
    imgsubhdr.IID2 = strtrim(fread(fid,80,'uint8=>char')');
    imgsubhdr.ISCLAS = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISCLSY = strtrim(fread(fid,2,'uint8=>char')');
    imgsubhdr.ISCODE = strtrim(fread(fid,11,'uint8=>char')');
    imgsubhdr.ISCTLH = strtrim(fread(fid,2,'uint8=>char')');
    imgsubhdr.ISREL = strtrim(fread(fid,20,'uint8=>char')');
    imgsubhdr.ISDCTP = strtrim(fread(fid,2,'uint8=>char')');
    imgsubhdr.ISDCDT = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ISDCXM = strtrim(fread(fid,4,'uint8=>char')');
    imgsubhdr.ISDG = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISDGDT = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ISCLTX = strtrim(fread(fid,43,'uint8=>char')');
    imgsubhdr.ISCATP = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISCAUT = strtrim(fread(fid,40,'uint8=>char')');
    imgsubhdr.ISCRSN = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISSRDT = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ISCTLN = strtrim(fread(fid,15,'uint8=>char')');
    imgsubhdr.ENCRYP = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISORCE = strtrim(fread(fid,42,'uint8=>char')');
    imgsubhdr.NROWS = str2double(fread(fid,8,'uint8=>char')');
    imgsubhdr.NCOLS = str2double(fread(fid,8,'uint8=>char')');
    imgsubhdr.PVTYPE = strtrim(fread(fid,3,'uint8=>char')');
    imgsubhdr.IREP = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ICAT = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ABPP = str2double(fread(fid,2,'uint8=>char')');
    imgsubhdr.PJUST = fread(fid,1,'uint8=>char')';

    imgsubhdr.ICORDS = strtrim(fread(fid,1,'uint8=>char')');
    if ~strcmp(imgsubhdr.ICORDS,'')
        imgsubhdr.IGEOLO = fread(fid,60,'uint8=>char')';
    end

    imgsubhdr.NICOM = str2double(fread(fid,1,'uint8=>char')');
    for lp = 1:imgsubhdr.NICOM
        fld1 = strcat('ICOM',num2str(lp));
        imgsubhdr.(fld1) = strtrim(fread(fid,80,'uint8=>char')');
    end

    imgsubhdr.IC = fread(fid,2,'uint8=>char')';
    if ~strcmp(imgsubhdr.IC,'NC') && ~strcmp(imgsubhdr.IC,'NM')
        imgsubhdr.COMRAT = strtrim(fread(fid,4,'uint8=>char')');
    end

    imgsubhdr.NBANDS = str2double(fread(fid,1,'uint8=>char')');
    if imgsubhdr.NBANDS == 0
        imgsubhdr.XBANDS = str2double(fread(fid,5,'uint8=>char')');
        numBands = imgsubhdr.XBANDS;
    else
        numBands = imgsubhdr.NBANDS;
    end

    for lp = 1:numBands
        fld1 = strcat('IREPBAND',num2str(lp));
        imgsubhdr.(fld1) = strtrim(fread(fid,2,'uint8=>char')');
        fld2 = strcat('ISUBCAT',num2str(lp));
        imgsubhdr.(fld2) = strtrim(fread(fid,6,'uint8=>char')');
        fld3 = strcat('IFC',num2str(lp));
        imgsubhdr.(fld3) = fread(fid,1,'uint8=>char')';
        fld4 = strcat('IMFLT',num2str(lp));
        imgsubhdr.(fld4) = strtrim(fread(fid,3,'uint8=>char')');

        fld5 = strcat('NLUTS',num2str(lp));
        imgsubhdr.(fld5) = str2double(fread(fid,1,'uint8=>char')');
        if imgsubhdr.(fld5) > 0
            fld6 = strcat('NELUT',num2str(lp));
            imgsubhdr.(fld6) = str2double(fread(fid,5,'uint8=>char')');
            for numlut = 1:imgsubhdr.(fld5)  % Number of bands
                fldN = strcat('LUTD',num2str(numlut));
                for ct = 1:imgsubhdr.(fld6)  % Entries per band
                    fld7 = strcat('LUTD',num2str(ct-1));
                    imgsubhdr.LUTD.(fldN).(fld7) = uint8(fread(fid,1,'uint8=>char')');
                end
            end
        end
    end

    imgsubhdr.ISYNC = fread(fid,1,'uint8=>char')';
    imgsubhdr.IMODE = fread(fid,1,'uint8=>char')';
    imgsubhdr.NBPR = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NBPC = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NPPBH = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NPPBV = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NBPP = str2double(fread(fid,2,'uint8=>char')');
    imgsubhdr.IDLVL = str2double(fread(fid,3,'uint8=>char')');
    imgsubhdr.IALVL = str2double(fread(fid,3,'uint8=>char')');
    imgsubhdr.ILOC = fread(fid,10,'uint8=>char')';
    imgsubhdr.IMAG = fread(fid,4,'uint8=>char')';
    if strcmp(imgsubhdr.IMAG(1),'/')
        imgsubhdr.IMAG = 1.0 / str2double(imgsubhdr.IMAG(2:4));
    else
        imgsubhdr.IMAG = str2double(imgsubhdr.IMAG);
    end

    imgsubhdr.UDIDL = str2double(fread(fid,5,'uint8=>char')');
    if imgsubhdr.UDIDL > 0
        imgsubhdr.UDOFL = str2double(fread(fid,3,'uint8=>char')');
    end

    imgsubhdr.IXSHDL = str2double(fread(fid,5,'uint8=>char')');
    if imgsubhdr.IXSHDL > 0
        imgsubhdr.IXSOFL = str2double(fread(fid,3,'uint8=>char')');
    end

%% VERSION 2.0

elseif strcmp(nitfver,'2.0')
    imgsubhdr.IM = fread(fid,2,'uint8=>char')';
    imgsubhdr.IID = strtrim(fread(fid,10,'uint8=>char')');
    imgsubhdr.IDATIM = fread(fid,14,'uint8=>char')';
    imgsubhdr.TGTID = strtrim(fread(fid,17,'uint8=>char')');
    imgsubhdr.ITITLE = strtrim(fread(fid,80,'uint8=>char')');
    imgsubhdr.ISCLAS = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISCODE = strtrim(fread(fid,40,'uint8=>char')');
    imgsubhdr.ISCTLH = strtrim(fread(fid,40,'uint8=>char')');
    imgsubhdr.ISREL = strtrim(fread(fid,40,'uint8=>char')');
    imgsubhdr.ISCAUT = strtrim(fread(fid,20,'uint8=>char')');
    imgsubhdr.ISCTLN = strtrim(fread(fid,20,'uint8=>char')');
    imgsubhdr.ISDWNG = strtrim(fread(fid,6,'uint8=>char')');

    if strcmp(imgsubhdr.ISDWNG,'999998')
        imgsubhdr.ISDEVT = strtrim(fread(fid,40,'uint8=>char')');
    end

    imgsubhdr.ENCRYP = fread(fid,1,'uint8=>char')';
    imgsubhdr.ISORCE = strtrim(fread(fid,42,'uint8=>char')');
    imgsubhdr.NROWS = str2double(fread(fid,8,'uint8=>char')');
    imgsubhdr.NCOLS = str2double(fread(fid,8,'uint8=>char')');
    imgsubhdr.PVTYPE = fread(fid,3,'uint8=>char')';
    imgsubhdr.IREP = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ICAT = strtrim(fread(fid,8,'uint8=>char')');
    imgsubhdr.ABPP = str2double(fread(fid,2,'uint8=>char')');
    imgsubhdr.PJUST = fread(fid,1,'uint8=>char')';

    imgsubhdr.ICORDS = strtrim(fread(fid,1,'uint8=>char')');
    if ~strcmp(imgsubhdr.ICORDS,'') && ~strcmp(imgsubhdr.ICORDS,'N')
        imgsubhdr.IGEOLO = fread(fid,60,'uint8=>char')';
    end

    imgsubhdr.NICOM = str2double(fread(fid,1,'uint8=>char')');
    for lp = 1:imgsubhdr.NICOM
        fld1 = strcat('ICOM',num2str(lp));
        imgsubhdr.ICOM.(fld1) = strtrim(fread(fid,80,'uint8=>char')');
    end

    imgsubhdr.IC = fread(fid,2,'uint8=>char')';
    if ~strcmp(imgsubhdr.IC,'NC') && ~strcmp(imgsubhdr.IC,'NM')
        imgsubhdr.COMRAT = strtrim(fread(fid,4,'uint8=>char')');
    end

    imgsubhdr.NBANDS = str2double(fread(fid,1,'uint8=>char')');
    for lp = 1:imgsubhdr.NBANDS
        fld1 = strcat('IREPBAND',num2str(lp));
        imgsubhdr.(fld1) = strtrim(fread(fid,2,'uint8=>char')');
        fld2 = strcat('ISUBCAT',num2str(lp));
        imgsubhdr.(fld2) = strtrim(fread(fid,6,'uint8=>char')');
        fld3 = strcat('IFC',num2str(lp));
        imgsubhdr.(fld3) = fread(fid,1,'uint8=>char')';
        fld4 = strcat('IMFLT',num2str(lp));
        imgsubhdr.(fld4) = strtrim(fread(fid,3,'uint8=>char')');

        fld5 = strcat('NLUTS',num2str(lp));
        imgsubhdr.(fld5) = str2double(fread(fid,1,'uint8=>char')');
        if imgsubhdr.(fld5) > 0
            fld6 = strcat('NELUT',num2str(lp));
            imgsubhdr.(fld6) = str2double(fread(fid,5,'uint8=>char')');
            for numlut = 1:imgsubhdr.(fld5)  % Number of bands
                fldN = strcat('LUTD',num2str(numlut));
                for ct = 1:imgsubhdr.(fld6)  % Entries per band
                    fld7 = strcat('LUTD',num2str(ct-1));
                    imgsubhdr.LUTD.(fldN).(fld7) = uint8(fread(fid,1,'uint8=>char')');
                end
            end
        end
    end

    imgsubhdr.ISYNC = fread(fid,1,'uint8=>char')';
    imgsubhdr.IMODE = fread(fid,1,'uint8=>char')';
    imgsubhdr.NBPR = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NBPC = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NPPBH = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NPPBV = str2double(fread(fid,4,'uint8=>char')');
    imgsubhdr.NBPP = str2double(fread(fid,2,'uint8=>char')');
    imgsubhdr.IDLVL = str2double(fread(fid,3,'uint8=>char')');
    imgsubhdr.IALVL = str2double(fread(fid,3,'uint8=>char')');
    imgsubhdr.ILOC = fread(fid,10,'uint8=>char')';
    imgsubhdr.IMAG = str2double(fread(fid,4,'uint8=>char')');

    imgsubhdr.UDIDL = str2double(fread(fid,5,'uint8=>char')');
    if imgsubhdr.UDIDL > 0
        imgsubhdr.UDOFL = str2double(fread(fid,3,'uint8=>char')');
        imgsubhdr.UDID = fread(fid,imgsubhdr.UDIDL-3,'uint8=>char')';
    end

    imgsubhdr.IXSHDL = str2double(fread(fid,5,'uint8=>char')');
    if imgsubhdr.IXSHDL > 0
        imgsubhdr.IXSOFL = str2double(fread(fid,3,'uint8=>char')');
    end
end

%% Read each extension
if imgsubhdr.IXSHDL > 0
    IXSHDLremainder = imgsubhdr.IXSHDL-3;
    while IXSHDLremainder > 0
        % Read TRE Tag and Length
        CETAG = strtrim(fread(fid,6,'uint8=>char')');
        CEL = str2double(fread(fid,5,'uint8=>char')');
        fseek(fid, -11, 'cof'); % Reset position indicator
        % Check to see if reader function exists for this TRE.
        if exist(strcat('read',CETAG,'.m'),'file');
            % Create CETAG structure from dynamic variable.
            tre  =  eval(strcat('read',CETAG,'(fid)'));
        else
            % Read undefined data
            tre = readUndefinedTRE(fid);
        end
        % Add dynamic field names to structure
        if ~isfield(imgsubhdr,CETAG)
            imgsubhdr.(CETAG) = tre;
        else
            imgsubhdr.(CETAG)(end+1,1) = tre;
        end
        % Calculate remaining number of subheader bytes
        IXSHDLremainder = IXSHDLremainder - (CEL + 11);
    end
end

%% Image data mask table
% Assure big-endian reads since it matters here
if any(imgsubhdr.IC=='M')
    imgsubhdr.IMDATOFF = fread(fid,1,'uint32',0,'b');
    imgsubhdr.BMRLNTH = fread(fid,1,'uint16',0,'b');
    imgsubhdr.TMRLNTH = fread(fid,1,'uint16',0,'b');
    imgsubhdr.TPXCDLNTH = fread(fid,1,'uint16',0,'b');
    if imgsubhdr.TPXCDLNTH>0
        bytelength=ceil(double(imgsubhdr.TPXCDLNTH)/8);
        switch imgsubhdr.PVTYPE
            case 'INT'
                imgsubhdr.TPXCD = fread(fid,1,['uint' num2str(bytelength*8)],0,'b');
            case 'SI'
                imgsubhdr.TPXCD = fread(fid,1,['int' num2str(bytelength*8)],0,'b');
            case 'R'
                imgsubhdr.TPXCD = fread(fid,1,['float' num2str(imgsubhdr.NBPP)],0,'b');
            case 'C'
                imgsubhdr.TPXCD = fread(fid,2,['float' num2str(imgsubhdr.NBPP)/2],0,'b');
        end
    end
    if imgsubhdr.BMRLNTH>0
        n=imgsubhdr.NBPR*imgsubhdr.NBPC;
        if imgsubhdr.IMODE=='S', m=imgsubhdr.NBANDS; else m=1; end;
        imgsubhdr.BMRnBNDm=zeros(n,m);
        for i=1:m
            for j=1:n
                imgsubhdr.BMRnBNDm(j,i) = fread(fid,1,'uint32',0,'b');
            end
        end
    end
    if imgsubhdr.TMRLNTH>0
        n=imgsubhdr.NBPR*imgsubhdr.NBPC;
        if imgsubhdr.IMODE=='S', m=imgsubhdr.NBANDS; else m=1; end;
        imgsubhdr.TMRnBNDm=zeros(n,m);
        for i=1:m
            for j=1:n
                imgsubhdr.TMRnBNDm(j,i) = fread(fid,1,'uint32',0,'b');
            end
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////