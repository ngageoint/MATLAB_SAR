function fileList = image_search(dirname, search_filter_handle)
% IMAGE_SEARCH Searches the specified directory for images satisfying the
% given filter.
%
% This function reads the metadata from each image in the specified
% directory and calls the provided function handle on each one.  If the
% (boolean) result of the function is true the file name and metadata are
% returned in the resulting list.  If it evaluates to false, the image is
% ignored.
%
% Inputs:
%       dirname:        The string name of the directory in which the
%                       images live.  All images will be opened, checked,
%                       and closed.
%
%       search_filter_handle:  The search filter to be applied to each
%                       image found.  This is a function pointer (handle)
%                       to a function that takes on SICD structure and
%                       returns a boolean value.  True means "keep the
%                       image", false means "ignore the image".
%
% Return:
%       IMAGE_SEARCH returns an array, each element being a structure
%       with elements two elements: 
%               filename: The full path to the image
%               meta:     The SICD metadata structure for the image
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fullFileNameList = {};
metadataList     = {};
wb_hand=waitbar(0,'Enumerating files to search');
allFilesToReview = recurdir(dirname);
for i = 1:length(allFilesToReview)
    waitbar((i-1)/length(allFilesToReview),wb_hand,'Searching for recognized formats');
    try
        % disp(allFilesToReview{i}); % For debugging
        filename = fullfile(dirname,allFilesToReview{i});
        if ~isempty(guess_ph_format(filename))
            newreader = open_ph_reader(filename);
            % Maybe should convert to SICD metadata here (meta2sicd_cphdx)
            % for consistent filtering?  For now, leave as CPHD.
        else
            newreader = open_reader(filename);
        end
        if iscell(newreader)
            % Multiple images in file.  Assume each is an independant image and
            % tack on its metadata
            for iii=1:length(newreader)
                metadataList{end+1}     = newreader{iii}.get_meta();
                if isfield(metadataList{end},'native')
                    metadataList{end}=rmfield(metadataList{end},'native'); % Helps clear memory
                end
                fullFileNameList{end+1} = allFilesToReview{i};
                newreader{iii}.close();
            end
        else
            metadataList{end+1}     = newreader.get_meta();
            if isfield(metadataList{end},'native')
                metadataList{end}=rmfield(metadataList{end},'native'); % Helps clear memory
            end
            fullFileNameList{end+1} = allFilesToReview{i};
            newreader.close();
        end
    catch % Nothing to do.  Not a recognized format.  Just go to next file
    end
    clear newreader; % Helps clear memory
end


% Build the output array...
% Check each image to see if it passes the user-supplied filter.  If it
% does we'll append the file name to the end of the returned list,
% otherwise we'll just ignore/drop it.
fileList = [];
for i = 1:length(metadataList)
    filename = fullFileNameList{i};
    sicdmeta = metadataList{i};
    try
        if search_filter_handle(sicdmeta)
            result.filename = filename;
            result.meta     = sicdmeta;
            fileList = [fileList result];
        end
    end
    waitbar(1,wb_hand,'Filtering results');
end
close(wb_hand);

end

%RECURDIR List all files found recursively under given directory
function [ filelist ] = recurdir( dir2search )

if(nargin<1||isempty(dir2search))
    d = dir;
    dir2search=[];
else
    d = dir(dir2search);
end

filelist={};
for j=1:length(d)
    if ~d(j).isdir
        filelist{end+1,1}=d(j).name; % Record this file
    elseif ~strcmp(d(j).name(1),{'.','..','.svn'})
        sublist=recurdir(fullfile(dir2search, d(j).name)); % Find all of the files under this subdirectory
        sublist=cellfun(@(x) [d(j).name filesep x],sublist,'UniformOutput',false); % Append current directory
        filelist=cat(1,filelist,sublist); % Append to all previous files
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////