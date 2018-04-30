function [ filenames ] = palsar2_files( filename, check_nonstandard_filenames )
%PALSAR2_FILES Given one file in an ALOS PALSAR 2 data bundle of files, find the rest
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('check_nonstandard_filenames','var')
    % Should we look at other files in the directory by content, rather
    % than filename?
    check_nonstandard_filenames = false;
end

% A standard ALOS PALSAR 2 distribution will have several files all in the
% same directory.

FILETYPE_PREFIX = {'VOL-', 'LED-', 'IMG-[H|V]{2}-', 'TRL-'};
FILETYPE = cellfun(@(x) x(1:3), FILETYPE_PREFIX, 'UniformOutput', false);
RECTYPES.VOL = [1, 192, 192, 18, 18, 360];
RECTYPES.LED = [1, 11, 192, 18, 18, 720];
RECTYPES.IMG = [1, 50, 192, 18, 18, 720];
RECTYPES.TRL = [1, 63, 192, 18, 18, 720];

[pathstr, fnameonly, ext] = fileparts(filename);
fnameonly = [fnameonly ext];
dir_info = dir(pathstr);
for i = 1:numel(FILETYPE)
    filenames.(FILETYPE{i}) = {}; % Might be multiple
    % First check filename to see if its as expected.
    expected_filename = regexprep(fnameonly, ['^' strjoin(FILETYPE_PREFIX,'|')], FILETYPE_PREFIX(i));
    dir_info_type = dir_info(~cellfun(@isempty, regexp({dir_info.name}, expected_filename, 'once')));
    for j = 1:numel(dir_info_type)
        if check_filetype(fullfile(pathstr, dir_info_type(j).name), FILETYPE{i})
            filenames.(FILETYPE{i}){end+1} = fullfile(pathstr, dir_info_type(j).name);
        end
    end
    if check_nonstandard_filenames && isempty(filenames.(FILETYPE{i}))
        % If no expected filenames, check all files in current directory
        % for matching file type.
        % Warning: Could break if multiple datasets in the same
        % directory since we can't tell whether files are from same
        % dataset or not.  We assume files are either named
        % appropriately or grouped in directories appropriately, but
        % not necessarily both.
        for j = 1:numel(dir_info)
            if ~dir_info(j).isdir && ... 
                    check_filetype(fullfile(pathstr, dir_info(j).name), FILETYPE{i})
                filenames.(FILETYPE{i}){end+1} = fullfile(pathstr, dir_info(j).name);
            end
        end
    end
    % Cleanup
    if isempty(filenames.(FILETYPE{i}))
        filenames = rmfield(filenames, FILETYPE{i});
    elseif ~strcmp(FILETYPE{i},'IMG') % No cell arrays needed for non-IMG types
        filenames.(FILETYPE{i}) = filenames.(FILETYPE{i}){1};
    end
end

    function [ bool ] = check_filetype( filename, type )
        % Read volume descriptor records
        fid=fopen(filename,'r','b');
        cur_types = [fread(fid,1,'uint32') ... % Record sequence number
            fread(fid,1,'uint8') ... % Record subtype code
            fread(fid,1,'uint8') ... % Record type
            fread(fid,1,'uint8') ... % Record subtype code
            fread(fid,1,'uint8') ... % Record subtype code
            fread(fid,1,'uint32')]; % Record length
        fclose(fid);

        bool = isequal(cur_types, RECTYPES.(type));
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////