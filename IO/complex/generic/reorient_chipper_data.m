function [ reoriented_data ] = reorient_chipper_data( symmetry, data_in)
%REORIENT_CHIPPER_ARGS Applies a transform to the chipper indexing
%
% This function transforms the actual data that a chipper function outputs
% based on the symmetry specified.
%
% A chipper function is a simplified function of a file reader.  It takes
% only three arguments: [minIndexInDimension1 maxIndexInDimension1],
% [minIndexInDimension2 maxIndexInDimension2], [subsampleInDimension1,
% subsampleInDimension2].  All other information such as the file handle
% and other auxilliary information must be contained within the chipper
% function.
%
% The SYMMETRY parameter is used to perform symmetry operations on the
% imported image data.  A combination of mirroring in the first dimension
% (m1), mirroring in the second dimension (m2) and transpose(tr) is used to
% accomplish the symmetry operations.  The 'symmetry' input parameter is
% defined as [ 1m m2 tr ], where each of the elements can be either 0 or 1
% to indicate whether to apply a transformation or not.
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

%% Transform data
reoriented_data=data_in;
if symmetry(1)
    reoriented_data=flipud(reoriented_data);
end
if symmetry(2)
    reoriented_data=fliplr(reoriented_data);
end
if symmetry(3)
    reoriented_data=permute(reoriented_data,[2 1 3:ndims(reoriented_data)]);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////