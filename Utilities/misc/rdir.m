function [ filelist ] = rdir(searchdir,extensions)
%RDIR Recursive directory search with optional file extension filters
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% INPUTS:
%   searchdir: required - Search Directory 
%   extensions: optional - cell array of extensions filters
%
% OUTPUTS:
%   filelist: cell array of path/file names
%
% VERSION:
%   1.0  20110803  Tim Cox
%       - initial version, liberally borrowed some code from "image_search"
%         function (KRAUSS/SCHWARTZKOPF)
%
% TODO:

filelist={};

%make sure there is a filesep of the end of the searchdir
if ~(searchdir(end) == filesep)
    searchdir = [searchdir filesep];
end

%if extensions variable is not passed in then just return all file types
extensionsflag = 1;
if ~exist('extensions','var') ||  strcmp(extensions{1},'.*') 
    extensions = {'.*'};
    extensionsflag = 0;
else
    %make sure each extension has a dot as the first character
    for i=1:length(extensions)
        if ~(extensions{i}(1) == '.')
            extensions{i} = ['.' extensions{i}];
        end
    end
end

d = dir(searchdir);

for j=1:length(d)
    if ~d(j).isdir
        [pathstr, name, ext] = fileparts(d(j).name);
        foo = strcmp(ext,extensions);
        if max(foo(:)) == 1 || extensionsflag == 0
            filelist{1,end+1}=strcat(searchdir,d(j).name); % Record this file
        end
    elseif ~strcmp(d(j).name(1),{'.','..','.svn'})
        sublist=rdir([searchdir d(j).name filesep],extensions); % Find all of the files under this subdirectory
        %filter results to only include extensions specified
        sublist2 = [];
        count = 0;
        for i=1:length(sublist)
            [pathstr, name, ext] = fileparts(sublist{i});
            foo = strcmp(ext,extensions);
            if max(foo(:)) == 1 || extensionsflag == 0
                count = count+1;
                sublist2{count} = strcat('',sublist{i});
            end
        end
        filelist=horzcat(filelist,sublist2); % Append to all previous files
    end
end    

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////