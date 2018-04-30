function [ sicd_meta ] = sicd_update_meta( sicd_meta, version_str )
%SICD_UPDATE_META Update a SICD metadata structure to 1.0
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
    elseif compare_version(vers_parts, [0 5]) % Version 0.5
        sicd_meta = sicd_update_meta_0_5(sicd_meta);
    elseif compare_version(vers_parts, [0 4]) % Version 0.4
        sicd_meta = sicd_update_meta_0_4(sicd_meta);
    else % Either older version or misformed version string
        warning('SICD_UPDATE_META:SICD_VERSION','Unrecognized SICD version.  Behavior may not be as expected.');
        sicd_meta = sicd_update_meta_0_4(sicd_meta); % Attempt to update what we can
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