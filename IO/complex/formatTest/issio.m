function boolout = issio(filename)
% Streaming input/output format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[ihdr, endian, data_offset]=read_sio_meta(filename);
boolout=~isempty(ihdr);
if boolout
    dirout=dir(filename);
    % ihdr is uint32.  Without double cast, file size computations over 4 Gig can overflow
    if dirout.bytes ~= (double(data_offset) + (double(ihdr(2)) * double(ihdr(3)) * double(ihdr(5))));
        warning('ISSIO:malformedSIO',...
            'File appears to be SIO, but header does not match file size.');
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////