function [ metadata ] = read_caspr_meta( filename )
%READ_CASPR_META Read metadata from CASPR header
%   metadata = read_caspr_meta( filename )
%
%   Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid=fopen(filename,'r','b');
if(fid<0)
    error('READ_CASPR_META:InvalidFile','Unable to open file.');
end

% Read in header data
current_subfield='';
reading_subfield=false;
while ~feof(fid)
    linetoparse=fgetl(fid);
    if isempty(linetoparse)||linetoparse(1)==' ' % Empty field
        continue;
    end
    if strncmp(linetoparse,';;;',3)
        if ~isempty(regexprep(linetoparse,'[^a-zA-Z0-9]','')) %this check accounts for poorly formed RPHD files
            current_subfield=genvarname(deblank(regexprep(linetoparse,'[^a-zA-Z0-9]','')));
            reading_subfield = true;
        else
            reading_subfield=~reading_subfield;
            if reading_subfield
                current_subfield='';
            else
                current_subfield=current_subfield(1:min(end,63)); % Max length for MATLAB field names
            end
        end
        continue;
    end
    [value,fieldname]=strtok(linetoparse,'"'); % Some values with spaces are surrounded by quotes
    if isempty(fieldname)
        [value,fieldname]=strtok(linetoparse); % If not using quotes
    else
        fieldname=deblank(fieldname(2:end));
    end
    % We hardcode one field.
    [first, last] = regexp(upper(fieldname),'^\s*CLASSIFICATION\s*:?\s*');
    if any(first) && ~isempty(value) && value(1)==';'
        metadata.classification = fieldname(last+1:end);
    end
    % The rest of the fields are handled generically.
    if ~isempty(fieldname)
        fieldname=regexprep(fieldname,'[^a-zA-Z0-9]','');
        if reading_subfield % Subsection heading
            current_subfield=[current_subfield genvarname(deblank(fieldname))];
        elseif ~isempty(current_subfield) % Actual field value
            valid_fieldname=genvarname(fieldname);
            doubval=str2double(value);
            if(isnan(doubval)) % If value is numeric, store as such
                metadata.(current_subfield).(valid_fieldname)=value;
            else
                metadata.(current_subfield).(valid_fieldname)=doubval;
            end
        end
    end
end
fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////