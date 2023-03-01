function [ readerobj ] = open_csm_reader( filename )
%OPEN_CSM_READER Intiates a reader object for Cosmo Skymed HDF5 file format.
%
% This reader works for spotlight and stripmap (HIMAGE and PINGPONG) modes,
% but not SCANSAR modes.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Open HDF5 ids
fid=H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
mission_id = deblank(get_hdf_attribute(fid,'Mission ID')');
switch mission_id
    case {'CSK','KMPS'}
        dataset_str = '/SBI';
    otherwise  % 'CSG'
        dataset_str = '/IMG';
end
meta=meta2sicd_csm(fid);
num_bands=length(meta);
for i=1:num_bands % "pingpong" mode has multiple polarizations
    dset_id(i)=H5D.open(fid,['/S0' num2str(i) dataset_str]);
    dspace_id(i)=H5D.get_space(dset_id(i));
end
bands_closed=false(num_bands,1); % Remember which HDF data ids are open

%% Setup reader parameters
% Get datasize
[num_dims datasize] = H5S.get_simple_extent_dims(dspace_id(1)); % All polarizations should be same size
% Compute symmetry
symmetry=[0 0 1]; % Default for CSM complex data
lineorder=get_hdf_attribute(fid,'Lines Order')';
columnorder=get_hdf_attribute(fid,'Columns Order')';
% Should be EARLY-LATE/NEAR-FAR for all CSM complex data
symmetry(1)=~strncmp(columnorder,'NEAR-FAR',8);
symmetry(2)=xor(strncmp(lineorder,'EARLY-LATE',10),meta{1}.SCPCOA.SideOfTrack=='R');

%% Define object methods
readerobj=cell(num_bands,1);
for i=1:num_bands
    readerobj{i} = chipfun2readerobj(@(varargin) chipper_csm(i,varargin{:}), datasize);
    readerobj{i}.get_meta=@() meta{i};
    readerobj{i}.close=@() close_csm(i);
end

    function out = chipper_csm(band_index, varargin)
        [dim1range,dim2range,subsample]=check_chipper_args(datasize, varargin{:});
        tempargs=reorient_chipper_args(symmetry, datasize([2 1]), dim1range, dim2range, subsample);
        [dim1range,dim2range,subsample]=deal(tempargs{:});
        offset=[dim2range(1) dim1range(1)]-1;
        subsample([1 2]) = subsample([2 1]);
        slabsize=ceil(([diff(dim2range) diff(dim1range)]+1)./subsample);
        if num_dims>2 % Complex data
            offset(3)=0; subsample(3)=1; slabsize(3)=datasize(3);
        end
        H5S.select_hyperslab(dspace_id(band_index),'H5S_SELECT_SET',...
            offset,subsample,slabsize,ones(size(offset)));
        memspace_id = H5S.create_simple(length(slabsize),slabsize,slabsize);
        temp_out=H5D.read(dset_id(band_index),'H5ML_DEFAULT',memspace_id,dspace_id(band_index),'H5P_DEFAULT');
        H5S.close(memspace_id);
        if num_dims>2 % Complex data
            temp_out=squeeze(complex(temp_out(1,:,:),temp_out(2,:,:)));
        end
        out=reorient_chipper_data(symmetry,temp_out);
    end

    function close_csm(band_index)
        H5S.close(dspace_id(band_index));
        H5D.close(dset_id(band_index));
        bands_closed(band_index)=true;
        if all(bands_closed)
            H5F.close(fid);
        end
    end

    function value = get_hdf_attribute(hid_t,attribute_name)
        attribute=H5A.open_name(hid_t,attribute_name);
        value = H5A.read(attribute,H5A.get_type(attribute));
        H5A.close(attribute);
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////