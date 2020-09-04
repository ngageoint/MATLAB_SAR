function [ tags ] = read_tiff_tags( filename )
%READ_TIFF_TAGS Reads the tags from a TIFF file
%
% Returns a cell array of tags, each element of which contains a structure
% with the tags for each image in the TIFF file.
%
% Similar in behavior to MATLAB's IMFINFO, but specific to TIFF.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

t = Tiff(filename);
tagnames = t.getTagNames;
tags = {};
while true
    S = struct();
    for j = 1:numel(tagnames)
        try
            S.(tagnames{j}) = getTag(t, tagnames{j});
        end
    end
    tags{end+1} = S;
    if ~t.lastDirectory, t.nextDirectory, else break, end
end
t.close;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////