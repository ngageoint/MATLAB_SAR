
function fileList = image_pair_search(dirname, search_filter_handle)
% IMAGE_PAIR_SEARCH Searches the specified directory for image pairs
% satisfying the given filter.
%
% This function reads the metadata from each image in the specified
% directory and calls the provided function handle on each pair (i.e.
% images taken pairwise).  If the (boolean) result of the function is
% true the file name and metadata are returned in the resulting list.  
% If it evaluates to false, the image is ignored.
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
%       IMAGE_PAIR_SEARCH returns an array, each element being a structure
%       with elements two elements: 
%               filename1: The full path to the first image
%               filename2: The full path to the second image
%               meta1:     The SICD metadata structure for the first image
%               meta2:     The SICD metadata structure for the second image
%
% Notes: 
%   - This routine makes a single pass through the list of files in the
%     directory.  If the comparison is "oriented" pairs may be missed.  
%     That is, if there are 3 images in the directory, pair will be 
%     checked as
%          1,2;  1,3;  2,3;
%     If for some reason the pair 2,1 was different than pair 1,2 the
%     comparison will not be made and the pair (e.g. 2,1) will not be found.
%
%   - Images are taken pair-wise, so if there are more than 2 that satisfy
%     the constraints multiple pairs will be returned.  Suppose, for
%     example, that all three images all satisfy the constraint.  A list
%     than contains the pairs 1,2; 1,3; 2,3 will be returned.
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fullFileNameList = {};
metadataList     = {};
count            = 0;
D                = dir(dirname);
for i = 1:length(D)
  filename = D(i).name;
  if isdir(filename)
    continue
  end
  
  fullFileName = [dirname '/' filename];
  try
    newreader = open_reader(fullFileName);
    if iscell(newreader)
      % Multiple images in file.  Assume each is an independant image and
      % tack on its metadata
      for iii=1:size(newreader)
        count = count + 1;
        metadataList{count}     = newreader{iii}.get_meta();
        fullFileNameList{count} = fullFileName;
      end
    else
      count = count + 1;
      metadataList{count}     = newreader.get_meta();
      fullFileNameList{count} = fullFileName;
    end
    newreader.close();
  catch e
  end
end


% Build the output array...
% Check each image to see if it passes the user-supplied filter.  If it
% does we'll append the file name to the end of the returned list,
% otherwise we'll just ignore/drop it.
fileList = [];
for i = 1:length(metadataList)
  filename1 = fullFileNameList{i};
  sicdmeta1 = metadataList{i};
  for j = (i+1):length(metadataList)
    filename2 = fullFileNameList{j};
    sicdmeta2 = metadataList{j};
    if search_filter_handle(sicdmeta1,sicdmeta2)
      result.filename1 = filename1;
      result.meta1     = sicdmeta1;
      result.filename2 = filename2;
      result.meta2     = sicdmeta2;
      fileList = [fileList result]; 
    end
  end
end


end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
