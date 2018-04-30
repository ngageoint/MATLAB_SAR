function meta = read_cphd_preamble(filename)
% READ_CPHD_PREAMBLE Read preamble metadata from Compensated Phase History
% Data (CPHD) file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(filename, 'r');

meta = struct();
nb_line = '';
subfieldname='';
while ~ strncmp(nb_line, 'EndPreamble', 11) % Read metadata until "EndPreamble"
    nb_line = fgetl(fid);
    if strncmp(nb_line, 'BegPreamble',11), continue; end;
    [tag, value]=strtok(nb_line); % Break into tag and value parts
    tag=strtrim(tag);  value=strtrim(value); % Remove whitespace
    if ~isempty(tag)
        if strncmp(tag,'Beg',3) % Beg*/End* structure
            % Record per-pulse metadata fieldnames
            if strncmp(tag,'BegNBVector',11)||strncmp(tag,'BegNBPulse',10)
                % Make a CPHDX style VectorParameters field
                vbfieldname = strtrim(strtok(fgetl(fid)));
                while ~strncmp(vbfieldname, 'End', 3)
                    if strncmpi(vbfieldname,'SRP',3)||strcmpi(vbfieldname(end-2:end),'Pos')
                        meta.VectorParameters.(vbfieldname) = 24;
                    elseif strfind(upper(vbfieldname),'NUMBER') == length(vbfieldname)-5
                        meta.VectorParameters.(vbfieldname) = 4;
                    else
                        meta.VectorParameters.(vbfieldname) = 8;
                    end
                    vbfieldname = strtrim(strtok(fgetl(fid)));
                end
            else % Skip BegWBVector/EndWBVector structure info
                while ~strncmp(fgetl(fid), 'End', 3), end;
            end
            subfieldname='';
        elseif isempty(value) % Structured parameter
            subfieldname=tag;
        elseif strncmp(tag,'DateTime',8) % Handle date/time specifically
            try
                meta.(tag) = datenum(value,'yyyymmddHHMMSS');
            end
        else % Everything else is handled generically
            doublevalue=str2double(value); % Is value string or numeric
            if(isfinite(doublevalue)) % If value is numeric, store as such
                value=doublevalue;
            else
                switch value % Store booleans as booleans
                    case 'Yes'
                        value=true;
                    case 'No'
                        value=false;
                end
            end
            if isempty(subfieldname)
                meta.(tag) = value;
            else
                meta.(subfieldname).(tag) = value;
            end
        end
    end
end

magic_number = fread(fid, 1, 'uint32=>uint32','b');
switch magic_number
    case hex2dec('4D6F6A6F') % Big endian
        meta.endian='b';
    case hex2dec('6F6A6F4D') % Little endian
        meta.endian='l';
    otherwise % Didn't find magic code
        error('Not CPHD file: magic code not found');
end

meta.pulseStart=ftell(fid); % Record where pulse data starts in the file
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////