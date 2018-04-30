function [ new_input_args ] = reorient_chipper_args( symmetry, datasize, varargin)
%REORIENT_CHIPPER_ARGS Applies a transform to the chipper indexing
%
% This function does nothing but transform the arguments into a chipper
% function based on the symmetry specified.
%
% A chipper function is a simplified function of a file reader.  It takes
% only three arguments: [minIndexInDimension1 maxIndexInDimension1],
% [minIndexInDimension2 maxIndexInDimension2], [subsampleInDimension1,
% subsampleInDimension2].  All other information such as the file handle
% and other auxilliary information must be contained within the chipper
% function.
%
% The SYMMETRY parameter is used to perform symmetry operations on the
% imported image data.  A combination of horizontal mirroring (hm), vertical
% mirroring (vm) and transpose(tr) is used to accomplish the symmetry
% operations.  The 'symmetry' input parameter is defined as [ hm vm tr ], where
% each of the elements can be either 0 or 1 to indicate whether to apply a
% transformation or not.
%
%   [ 0 0 0 ] - default setting; successive pixels on disk are interpreted
%               to fill the first dimension in a MATLAB array (vertical, if
%               viewed in imagesc/imshow without any reorienting)
%   [ 1 0 1 ] -  90 degrees CCW from [ 0 0 0 ]
%   [ 0 1 1 ] -  90 degrees  CW from [ 0 0 0 ]
%   [ 1 1 0 ] - 180 degrees     from [ 0 0 0 ]
%
%   [ 0 0 1 ] - transpose of [ 0 0 0 ];
%   [ 0 1 0 ] -  90 degrees CCW from [ 0 0 1 ]
%   [ 1 0 0 ] -  90 degrees  CW from [ 0 0 1 ]
%   [ 1 1 1 ] - 180 degrees     from [ 0 0 1 ]
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if symmetry(3)
    new_input_args={};
    if(nargin>=3), new_input_args{2}=varargin{1}; end;
    if(nargin>=4), new_input_args{1}=varargin{2}; end
    if(nargin>=5), new_input_args{3}=varargin{3}(end:-1:1); end;
else
    new_input_args=varargin;
end
if symmetry(1)
    if(nargin>=3)
        new_input_args{1}=double(datasize(1))-new_input_args{1}(end:-1:1)+1;
        % The double around datasize allows this to work on empty arguments.
    end
end
if symmetry(2)
    if(nargin>=4)
        new_input_args{2}=double(datasize(2))-new_input_args{2}(end:-1:1)+1;
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////