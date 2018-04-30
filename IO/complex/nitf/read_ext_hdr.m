function extHdr = read_ext_hdr(fid,XHDL)
% READ_EXT_HDR Parse TREs from an extended header
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

extHdr=struct();
%% Read each TRE in extended header
if XHDL > 0
    XHDLremainder = XHDL-3;
    while XHDLremainder > 0
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
        if ~isfield(extHdr,CETAG)
            extHdr.(CETAG) = tre;
        else
            extHdr.(CETAG)(end+1,1) = tre;
        end
        % Calculate remaining number of subheader bytes
        XHDLremainder = XHDLremainder - (CEL + 11);
    end
end


end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////