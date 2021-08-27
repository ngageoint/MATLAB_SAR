function [ cphd_meta ] = cphd_update_meta( cphd_meta, version_str )
%CPHD_UPDATE_META Update a CPHD metadata structure to latest
%
% Basically a cut-and-paste from sicd_update_meta.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Update metadata structure to version 1.0 if necessary
vers_parts = str2double(regexp(version_str,'\.','split'));
if all(isfinite(vers_parts)) % Valid version format
    if compare_version(vers_parts, 1) % Version 1.0
        % do nothing
    elseif compare_version(vers_parts, [0 3]) % Version 0.3
        cphd_meta = cphd_update_meta_0_3(cphd_meta);
    else % Either older version or misformed version string
        error('CPHD_UPDATE_META:CPHD_VERSION','Unrecognized CPHD version.');
    end
end

end

% Is vers_ref is new than or equal to vers_comapre
function bool_out = compare_version(vers_ref, vers_compare)
    if numel(vers_ref)>numel(vers_compare)
        vers_compare(numel(vers_ref))=0;
    elseif numel(vers_compare)>numel(vers_ref)
        vers_ref(numel(vers_compare))=0;
    end
    tenpowers = 10.^(0:-1:-(numel(vers_ref)-1)).';
    bool_out = (sign(vers_compare - vers_ref) * tenpowers) <= 0;
end
    
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////