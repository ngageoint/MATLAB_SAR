function write_nitf_security_tags(obj)
%WRITE_NITF_SECURITY_TAGS Writes the NITF security tags
%
% It's not expected that the typical user would call this function, rather
% it is called as part of a bigger NITF/SICD writer.
%
% Inputs:
%       obj:       This is the 'hidden' object-specific handle.  It is
%                  similar to 'this' in Java.
%
% References:
%       MIL-STD-2500C, Department of Defense Interface Standard 
%                      National Imagery Transmission Format version 2.1
%                      01 May 2006
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Always writing version 2.1 NITF

% Security is limited in the SICD XML structure.  We fill in what we can.
if isfield(obj.sicdmeta,'CollectionInfo')&&isfield(obj.sicdmeta.CollectionInfo,'Classification')
    classification = obj.sicdmeta.CollectionInfo.Classification;
else
    classification = ' ';
end
code_index = regexp(classification,'/[^/]');
if isempty(code_index)
    code = '';
else
    code = classification((code_index(1)+1):end);
end
fwriten(obj.FID, classification(1), 1); % CLAS
fwriten(obj.FID, 'US',      2);  % CLSY
fwriten(obj.FID, code,      11); % CODE
fwriten(obj.FID, '',        2);  % CTLH
fwriten(obj.FID, '',        20); % REL
fwriten(obj.FID, '',        2);  % DCTP
fwriten(obj.FID, '',        8);  % DCDT
fwriten(obj.FID, '',        4);  % DCXM
fwriten(obj.FID, '',        1);  % DG
fwriten(obj.FID, '',        8);  % DGDT
fwriten(obj.FID, ''       , 43); % CLTX
fwriten(obj.FID, '',        1);  % CATP
fwriten(obj.FID, '',        40); % CAUT
fwriten(obj.FID, '',        1);  % CRSN
fwriten(obj.FID, '',        8);  % SRDT
fwriten(obj.FID, '',        15); % CTLN

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

