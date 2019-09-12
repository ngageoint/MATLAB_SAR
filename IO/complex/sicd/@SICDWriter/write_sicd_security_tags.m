function write_sicd_security_tags(obj)
%WRITE_SICD_SECURITY_TAGS Writes the NITF security tags
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

% Security is limited in the SICD XML structure.  Without explicit tags
% passed, we make a guess from XML and fill in what we can.
if isfield(obj.sicdmeta,'CollectionInfo') && ...
        isfield(obj.sicdmeta.CollectionInfo,'Classification')&& ...
        ~isempty(obj.sicdmeta.CollectionInfo.Classification)
    classification = obj.sicdmeta.CollectionInfo.Classification;
else
    classification = ' ';
end
codes = strsplit(classification,'/');  % Cell array
codes = strjoin(cellfun(@(x) x(1:min(2,end)), codes(2:end), 'UniformOutput', false));

% Tag name, field length, default value
tag_info = { ...
    'CLAS', 1, classification(1); ...
    'CLSY', 2, 'US'; ...
    'CODE', 11, codes; ...
    'CTLH', 2, ''; ...
    'REL', 20, ''; ...
    'DCTP', 2, ''; ...
    'DCDT', 8, ''; ...
    'DCXM', 4, ''; ...
    'DG', 1, ''; ...
    'DGDT', 8, ''; ...
    'CLTX', 43, ''; ...
    'CATP', 1, ''; ...
    'CAUT', 40, ''; ...
    'CRSN', 1, ''; ...
    'SRDT', 8, ''; ...
    'CTLN', 15, ''};

if isfield(obj.sicdmeta,'NITF') && isfield(obj.sicdmeta.NITF,'security')
    sec_tags = obj.sicdmeta.NITF.security;
else
    sec_tags = struct();
end

for i=1:size(tag_info,1)
    if isfield(sec_tags, tag_info{i,1})
        fwriten(obj.FID, sec_tags.(tag_info{i,1}), tag_info{i,2}); 
    else
        fwriten(obj.FID, tag_info{i,3}, tag_info{i,2}); 
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////