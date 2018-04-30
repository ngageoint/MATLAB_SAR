function txtsubhdr = read_nitf_txtsubhdr(fid,nitfver)
%READ_NITF_TXTSUBHDR Create structure for required text subheader metadata
%fields.

% Written by Wade Schwartzkopf, NGA/IDT
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
    txtsubhdr.TE = fread(fid,2,'uint8=>char')';
    txtsubhdr.TEXTID = strtrim(fread(fid,7,'uint8=>char')');
    txtsubhdr.TXTALVL = strtrim(fread(fid,3,'uint8=>char')');
    txtsubhdr.TXTDT = fread(fid,14,'uint8=>char')';
    txtsubhdr.TXTITL = strtrim(fread(fid,80,'uint8=>char')');
    txtsubhdr.TSCLAS = fread(fid,1,'uint8=>char')';
    txtsubhdr.TSCLSY = strtrim(fread(fid,2,'uint8=>char')');
    txtsubhdr.TSCODE = strtrim(fread(fid,11,'uint8=>char')');
    txtsubhdr.TSCTLH = strtrim(fread(fid,2,'uint8=>char')');
    txtsubhdr.TSREL = strtrim(fread(fid,20,'uint8=>char')');
    txtsubhdr.TSDCTP = strtrim(fread(fid,2,'uint8=>char')');
    txtsubhdr.TSDCDT = strtrim(fread(fid,8,'uint8=>char')');
    txtsubhdr.TSDCXM = strtrim(fread(fid,4,'uint8=>char')');
    txtsubhdr.TSDG = fread(fid,1,'uint8=>char')';
    txtsubhdr.TSDGDT = strtrim(fread(fid,8,'uint8=>char')');
    txtsubhdr.TSCLTX = strtrim(fread(fid,43,'uint8=>char')');
    txtsubhdr.TSCATP = fread(fid,1,'uint8=>char')';
    txtsubhdr.TSCAUT = strtrim(fread(fid,40,'uint8=>char')');
    txtsubhdr.TSCRSN = fread(fid,1,'uint8=>char')';
    txtsubhdr.TSSRDT = strtrim(fread(fid,8,'uint8=>char')');
    txtsubhdr.TSCTLN = strtrim(fread(fid,15,'uint8=>char')');
    txtsubhdr.ENCRYP = fread(fid,1,'uint8=>char')';
    txtsubhdr.TXTFMT = strtrim(fread(fid,3,'uint8=>char')');
    txtsubhdr.TXSHDL = str2double(fread(fid,5,'uint8=>char')');
    if txtsubhdr.TXSHDL > 0
        txtsubhdr.TXSOFL = str2double(fread(fid,3,'uint8=>char')');
    end
%% VERSION 2.0
elseif strcmp(nitfver,'2.0')
    txtsubhdr.TXSHDL = 0; % Don't bother to read version 2.0
end

%% Read each extension
if txtsubhdr.TXSHDL > 0
    TXSHDLremainder = txtsubhdr.TXSHDL-3;
    while TXSHDLremainder > 0
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
        if ~isfield(txtsubhdr,CETAG)
            txtsubhdr.(CETAG) = tre;
        else
            txtsubhdr.(CETAG)(end+1,1) = tre;
        end
        % Calculate remaining number of subheader bytes
        TXSHDLremainder = TXSHDLremainder - (CEL + 11);
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////