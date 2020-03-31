function [ readerobj ] = open_nisar_reader( filename )
%OPEN_NISAR_READER Intiates a reader object for NISAR SLC HDF5 file format.
%
% Written by: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Open HDF5 ids
fid=H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
meta=meta2sicd_nisar(fid);

% Query for frequencies and polarization (each of which would create its
% own SICD).  Replicates the order used in meta2sicd_nisar.
freqs_id = H5D.open(fid,'/science/LSAR/identification/listOfFrequencies');
freqs = H5D.read(freqs_id);
freqs = freqs(1,:);
H5D.close(freqs_id);
for i=1:numel(freqs)
    pol_id = H5D.open(fid,['/science/LSAR/SLC/swaths/frequency' freqs(i) '/listOfPolarizations']);
    pols{i} = H5D.read(pol_id)';
    H5D.close(pol_id);
end

% Open HDF objects for reading complex pixels
meta_i=0;
for i=1:numel(freqs)
    for j=1:size(pols{i},1)
        meta_i = meta_i + 1;
        dset_id(meta_i)=H5D.open(fid,['/science/LSAR/SLC/swaths/frequency' freqs(i) '/' pols{i}(j,:)]);
        dspace_id(meta_i)=H5D.get_space(dset_id(meta_i));
        [~, datasize{meta_i}] = H5S.get_simple_extent_dims(dspace_id(meta_i));
    end
end
bands_closed=false(numel(meta),1); % Remember which HDF data ids are open

%% Setup reader parameters
% for i=1:8, [~,b] = H5S.get_simple_extent_dims(dspace_id(i)); disp(b), end
% Compute symmetry
% Plans are for NISAR to be left-looking
symmetry=[0 ~(meta{1}.SCPCOA.SideOfTrack=='R') 1];

%% Define object methods
readerobj=cell(numel(meta),1);
for i=1:numel(meta)
    readerobj{i} = chipfun2readerobj(@(varargin) chipper_nisar(i,varargin{:}), datasize{i});
    readerobj{i}.get_meta=@() meta{i};
    readerobj{i}.close=@() close_hdf(i);
end

    function out = chipper_nisar(band_index, varargin)
        [dim1range,dim2range,subsample]=check_chipper_args(datasize{band_index}, varargin{:});
        tempargs=reorient_chipper_args(symmetry, datasize{band_index}([2 1]), dim1range, dim2range, subsample);
        [dim1range,dim2range,subsample]=deal(tempargs{:});
        offset=[dim2range(1) dim1range(1)]-1;
        subsample([1 2]) = subsample([2 1]);
        slabsize=ceil(([diff(dim2range) diff(dim1range)]+1)./subsample);
        H5S.select_hyperslab(dspace_id(band_index),'H5S_SELECT_SET',...
            offset,subsample,slabsize,ones(size(offset)));
        memspace_id = H5S.create_simple(length(slabsize),slabsize,slabsize);
        temp_out=H5D.read(dset_id(band_index),'H5ML_DEFAULT',memspace_id,dspace_id(band_index),'H5P_DEFAULT');
        H5S.close(memspace_id);
        out=reorient_chipper_data(symmetry, temp_out.r + (1i*temp_out.i));
    end

    function close_hdf(band_index)
        H5S.close(dspace_id(band_index));
        H5D.close(dset_id(band_index));
        bands_closed(band_index)=true;
        if all(bands_closed)
            H5F.close(fid);
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////