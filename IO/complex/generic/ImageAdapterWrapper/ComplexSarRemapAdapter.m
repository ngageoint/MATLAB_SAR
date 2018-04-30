classdef ComplexSarRemapAdapter < ImageAdapter
% COMPLEXSARREMAPADAPTER Wraps the open_reader framework from the MATLAB
% SAR Toolbox in the MATLAB ImageAdapter class so that R-sets can be
% generated from complex SAR data.
%
% Currently a global 8-bit remap is applied.  This allows data to be viewed
% in imtool, but this means we lose dynamic range when viewing the image,
% and it prevents this class from being used for processing using MATLAB's
% blockproc.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

    properties(GetAccess = public, SetAccess = private)
        Filename
        Metadata
    end
    
    properties(GetAccess = private, SetAccess = private)
        reader_object;
        sample_mean;
    end
    
    methods
        function obj = ComplexSarRemapAdapter(filename)
            % Open file
            if ischar(filename) % Input can be single string or cell array of strings
                filename = {filename}; % Just treat both types as cell array
            end
            obj.reader_object = {};
            for i = 1:length(filename)
                obj.reader_object = [obj.reader_object open_reader(filename{i})];
            end
            [readerobj_indices,obj.reader_object]=group_reader_objs_by_pol(obj.reader_object);
            if numel(obj.reader_object)>1
                error('COMPLEXSARADAPTER:MULTIIMAGEFILE','ComplexSarAdapter does not currently handle multi-image datasets.');
            else
                obj.reader_object = obj.reader_object{1}; 
            end

            % Set properties
            obj.Filename = filename;
            obj.Metadata = obj.reader_object.get_meta();
            obj.ImageSize = double([obj.Metadata.ImageData.NumRows, ...
                obj.Metadata.ImageData.NumCols]);
            
            % Estimate data mean
            samplesize=[1000 1000]; % Exract 1000x1000 array of samples to estimate mean
            subsample=ceil(double([obj.Metadata.ImageData.NumCols obj.Metadata.ImageData.NumRows])./samplesize);
            data = abs(single(obj.reader_object.read_chip([1 obj.Metadata.ImageData.NumCols],...
                [1 obj.Metadata.ImageData.NumRows],subsample)));
            data = pol_decomp(data,obj.Metadata);
            for i = 1:size(data,3)
                single_band = data(:,:,i);
                obj.sample_mean(i)=mean(single_band(isfinite(single_band)));
            end
        end
        
        function data = readRegion(obj, region_start, region_size)
            range_1 = region_start(1) + [0, region_size(1) - 1];
            range_2 = region_start(2) + [0, region_size(2) - 1];
            data = obj.reader_object.read_chip(range_2, range_1);
            data = permute(data,[2 1 3:ndims(data)]); % Reorient to what imtool expects
            data = pol_decomp(data,obj.Metadata); % Apply polarimetric decomposition if necessary
            data = do_remap(data,obj.sample_mean); % Apply 8-bit remap
        end
        
        function close(obj)
            obj.reader_object.close();
        end
    end
end

% Temporarily cut-and-pasted from hg_mitm_viewer makeDisplayable
% Ugly hack to account for polarimetric and color data. In theory, if this
% function becomes useful, we should tie these two functions together.
function out = pol_decomp(in, metadata)
    out = in;

    if isinteger(out) % Many remap functions won't work on complex ints
        out=single(out);
    end

    switch size(out,3)
        case 1 % Single band image; nothing to do
        case 2 % Dual-pol (quasi-Pauli)
           co_index=[find(strcmpi('H:H',metadata.ImageFormation.TxRcvPolarizationProc))...
              find(strcmpi('V:V',metadata.ImageFormation.TxRcvPolarizationProc))];
           if length(co_index)>1 % Catch HH/VV
              cross_index=co_index(2);
              co_index=co_index(1);
           else
              cross_index=[find(strcmpi('H:V',metadata.ImageFormation.TxRcvPolarizationProc))...
                find(strcmpi('V:H',metadata.ImageFormation.TxRcvPolarizationProc))];
           end
           out=cat(3,abs(out(:,:,co_index)),...
              abs(out(:,:,cross_index)),...
              abs(out(:,:,co_index)));
        case 3 % RGB image; nothing to do
        case 4 % Quad-pol: default to Pauli decomposition
            HH_index=find(strcmpi('H:H',metadata.ImageFormation.TxRcvPolarizationProc));
            HV_index=find(strcmpi('H:V',metadata.ImageFormation.TxRcvPolarizationProc));
            VH_index=find(strcmpi('V:H',metadata.ImageFormation.TxRcvPolarizationProc));
            VV_index=find(strcmpi('V:V',metadata.ImageFormation.TxRcvPolarizationProc));
            out=cat(3,abs(out(:,:,HH_index)-out(:,:,VV_index)),...
                abs(out(:,:,HV_index)+out(:,:,VH_index))/2,...
                abs(out(:,:,HH_index)+out(:,:,VV_index)));
    end
end

function out = do_remap(in, sample_mean)
    % Apply remap per band
    for i=1:size(in,3)
        out(:,:,i)=uint8(amplitudetodensity(in(:,:,i),30,40,sample_mean(i)));
    end
end


% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////