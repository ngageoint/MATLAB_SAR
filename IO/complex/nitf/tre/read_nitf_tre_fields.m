function [ meta ] = read_nitf_tre_fields( fid, FormatSpec )
%READ_NITF_TRE_FIELDS Read a set of TRE fields described by FormatSpec
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

for i=1:size(FormatSpec,1)
    meta.(FormatSpec{i,1}) = fread(fid, FormatSpec{i,2}, 'uint8=>char')';
    if FormatSpec{i,3}
        meta.(FormatSpec{i,1})=str2double(meta.(FormatSpec{i,1}));
    end
end    

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////