function [ dim1range, dim2range, subsample ] =...
    check_chipper_args( datasize, dim1range, dim2range, subsample, varargin)
%CHECK_CHIPPER_ARGS Check arguments given to a chipper function against
%datasize for validity
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if((nargin<4)||isempty(subsample)), subsample=[1 1]; end;
if((nargin<3)||isempty(dim2range)), dim2range=[1 datasize(2)]; end;
if((nargin<2)||isempty(dim1range)), dim1range=[1 datasize(1)]; end;
if(dim2range(2)>datasize(2))||any(dim2range<=0)||(dim2range(1)>dim2range(2))||...
   (dim1range(2)>datasize(1))||any(dim1range<=0)||(dim1range(1)>dim1range(2))      
    error('CHECK_CHIPPER_ARGS:InvalidArg','Invalid subimage index range.');
end
if((nargin<1)||any(datasize<=0))
    error('CHECK_CHIPPER_ARGS:InvalidArg','Invalid datasize.');
end

% Cast all values to double.  uint32 values can result in overflows for
% large files (>4 Gig)
dim1range=double(dim1range);
dim2range=double(dim2range);
subsample=double(subsample);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////