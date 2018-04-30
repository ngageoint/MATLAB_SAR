function status = fseek( fid, offset, origin )
%FSEEK Set file position indicator. 
% STATUS = FSEEK(FID, OFFSET, ORIGIN) mimics the behavior of MATLAB's FSEEK
% function, but allows for seeking past the end of a file (just appends
% zeros to the end of the file until it reaches desired location).
%
%   OFFSET values are interpreted as follows:
%       >= 0   Move position indicator OFFSET bytes after ORIGIN.
%       < 0    Move position indicator OFFSET bytes before ORIGIN.
%
%   ORIGIN values are interpreted as follows:
%       'bof' or -1   Beginning of file
%       'cof' or  0   Current position in file
%       'eof' or  1   End of file
%
%   STATUS is 0 on success and -1 on failure.  If an error occurs, use
%   FERROR to get more information.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

status = fseek(fid, offset, origin); % Try MATLAB's fseek first
if status < 0 % Didn't work
    [message, errnum] = ferror(fid);
    if errnum == -27 % "Offset is bad - after end-of-file or last character written."
        % This is not an error in our case.  We will do special handling to
        % allow this fseek to skip to position past current end of file.
        switch origin
            case {-1, 'bof'}
                offset = offset - ftell(fid); % Convert offset to be from current position, not beginning of file
            case {0, 'cof'} % No change required
            case {1, 'eof'}
                fseek(fid,0,'eof');
        end
        if offset>0
            % Use fwrite skip parameter to append to end of file an
            % arbitrary amount
            status = -(fwrite(fid, uint8(0), 'uint8', offset-1) ~= 1);
        end
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////