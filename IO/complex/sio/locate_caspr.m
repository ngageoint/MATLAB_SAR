function [ caspr_filename ] = locate_caspr( sio_filename )
%LOCATE_CASPR Locate CASPR metadata file that might be associated with an
%SIO file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[path, name]=fileparts(sio_filename);
possible_filenames={fullfile(path,[name '.hydra']),...
    fullfile(path,[name '.hdr']),...
    fullfile(path,'..','RPHDHeader.out'),...
    fullfile(path,'RPHDHeader.out'),...
    fullfile(path,[name '_RPHDHeader.out'])};
caspr_filename='';
for i=1:length(possible_filenames)
    if exist(possible_filenames{i},'file')
        caspr_filename=possible_filenames{i};
        break;
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////