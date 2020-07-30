function [ readerobj ] = open_iceye_reader( filename )
%OPEN_CSM_READER Intiates a reader object for ICEYE SLC HDF5 file format.
%
% Written by: Tim Cox, NRL.  Based on the CSM/NISAR reader by Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Open HDF5
fid=H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
meta=meta2sicd_iceye(fid);

% Open HDF objects for reading real and imaginary pixels
dset_id(1) = H5D.open(fid,'s_i');
dspace_id(1) = H5D.get_space(dset_id(1));
[~, datasize{1}] = H5S.get_simple_extent_dims(dspace_id(1));
dset_id(2) = H5D.open(fid,'s_q');
dspace_id(2) = H5D.get_space(dset_id(2));
[~, datasize{2}] = H5S.get_simple_extent_dims(dspace_id(2));

if strcmpi(h5read(filename,'/look_side'),'left')
    symmetry=[0 1 1]; 
else
    symmetry=[0 0 1];
end

%% Define object methods
readerobj = chipfun2readerobj(@chipper_iceye, datasize{1}); %real/imag are the same size
readerobj.get_meta=@() meta;
readerobj.close=@() close_hdf;

    function out = chipper_iceye(varargin)
        [dim1range,dim2range,subsample]=check_chipper_args(datasize{1}, varargin{:});
        tempargs=reorient_chipper_args(symmetry, datasize{1}([2 1]), dim1range, dim2range, subsample);
        [dim1range,dim2range,subsample]=deal(tempargs{:});
        offset=[dim2range(1) dim1range(1)]-1;
        subsample([1 2]) = subsample([2 1]);
        slabsize=ceil(([diff(dim2range) diff(dim1range)]+1)./subsample);
        H5S.select_hyperslab(dspace_id(1),'H5S_SELECT_SET',...
            offset,subsample,slabsize,ones(size(offset)));
        memspacei_id = H5S.create_simple(length(slabsize),slabsize,slabsize);
        temp_out_real=H5D.read(dset_id(1),'H5ML_DEFAULT',memspacei_id,dspace_id(1),'H5P_DEFAULT');
        H5S.select_hyperslab(dspace_id(2),'H5S_SELECT_SET',...
            offset,subsample,slabsize,ones(size(offset)));
        memspaceq_id = H5S.create_simple(length(slabsize),slabsize,slabsize);
        temp_out_imag=H5D.read(dset_id(2),'H5ML_DEFAULT',memspaceq_id,dspace_id(2),'H5P_DEFAULT');
        
        H5S.close(memspacei_id);
        H5S.close(memspaceq_id);
        out=reorient_chipper_data(symmetry, complex(temp_out_real,temp_out_imag));
    end

    function close_hdf()
        H5S.close(dspace_id(1));
        H5D.close(dset_id(1));
        H5S.close(dspace_id(2));
        H5D.close(dset_id(2));
        H5F.close(fid);
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////