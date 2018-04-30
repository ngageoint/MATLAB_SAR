function [ readerobj ] = open_tiff_reader_noxml(filename, symmetry)
%OPEN_TIFF_READER Intiates a reader object for TIFF imagery
%
% Does NOT search for associated XML metadata files.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin<2
    symmetry=[0 0 0];
end

%% Setup local variables
meta.native.tiff=imfinfo(filename);
meta.ImageData.NumCols=uint32(meta.native.tiff(1).Height);
meta.ImageData.NumRows=uint32(meta.native.tiff(1).Width);
datasize=[meta.ImageData.NumCols meta.ImageData.NumRows];

%% Define object methods
readerobj = chipfun2readerobj(@chipper, datasize);
readerobj.get_meta=@() meta;

    function out = chipper(varargin)
        newargs=reorient_chipper_args(symmetry, datasize, varargin{:});
        [dim1range,dim2range,subsample]=check_chipper_args(datasize, newargs{:});
        out = imread(filename,'tiff', 'PixelRegion', ...
            {[dim1range(1) subsample(1) dim1range(2)],...
            [dim2range(1) subsample(2) dim2range(2)]});
        if(size(out,3)==2) % Complex image
            out=complex(out(:,:,1),out(:,:,2));
        end
        out=reorient_chipper_data(symmetry, out);
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////