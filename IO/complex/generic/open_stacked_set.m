function [ readerobj ] = open_stacked_set(varargin)
%OPEN_STACKED_SET Opens a set of reader objects and concatenates them to
%output a multiband (3-D) array
%
% readerobj = open_stacked_set('filename1','filename2',...)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

for i=1:length(varargin)
    ro{i}=open_reader(varargin{i});
    meta{i}=ro{i}.get_meta();
    if((meta{end}.ImageData.NumRows~=meta{1}.ImageData.NumRows)||...
            (meta{end}.ImageData.NumCols~=meta{1}.ImageData.NumCols))
        error('open_stacked_set:nonmatching_sizes','Sizes of all images must be the same.');
    end
end
readerobj=stack_readers(ro);
readerobj.get_meta=@() meta{1}; % Not sure how to merge metadata.  For now, just pick one.

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////