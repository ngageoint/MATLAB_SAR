function [ sicd_meta_current ] = sicd_update_meta_0_4( sicd_meta_0_4 )
%SICD_UPDATE_META_0_4 Update a SICD metadata structure from version 0.4 to
%current version (whatever that may be)
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

sicd_meta_0_5 = sicd_meta_0_4;
% Update WgtType format
for i = {'Row','Col'}
    if isfield(sicd_meta_0_4, 'Grid') && isfield(sicd_meta_0_4.Grid, i{1}) && ...
            isfield(sicd_meta_0_4.Grid.(i{1}),'WgtType') && ...
            ischar(sicd_meta_0_4.Grid.(i{1}).WgtType)
        wgt_name = sicd_meta_0_4.Grid.(i{1}).WgtType;
        [wgt_name, parameters] = strtok(wgt_name);
        sicd_meta_0_5.Grid.(i{1}).WgtType = struct('WindowName', wgt_name); % Reset from string
        if ~isempty(parameters)
            % Assumes only one parameter
            [parameter_name, parameter_value] = strtok(parameters,'= ');
            if numel(parameter_value)>1
                sicd_meta_0_5.Grid.(i{1}).WgtType.Parameter.name = parameter_name;
                sicd_meta_0_5.Grid.(i{1}).WgtType.Parameter.value = parameter_value(2:end);
            end
        end
    end
end
% We are now updated to version 0.5.  Now do rest of updates.
sicd_meta_current = sicd_update_meta_0_5(sicd_meta_0_5);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////