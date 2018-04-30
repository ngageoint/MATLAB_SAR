function [ chipper_function, close_function ] = generic_chipper( filename, datasize,...
    datatype, complextype, data_offset, endian, symmetry, bands_ip, blocksize)
%GENERIC_CHIPPER Summary of this function goes here
%
% A chipper function is a simplified function of a file reader.  It takes
% only three arguments: [minIndexInDimension1 maxIndexInDimension1],
% [minIndexInDimension2 maxIndexInDimension2], [subsampleInDimension1,
% subsampleInDimension2].  All other information such as the file handle
% and other auxilliary information must be contained within the chipper
% function.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup chipper
% Can we use memory mapping or fread/fwrite?
% If we are reading and writing from/to files in order, not much advantage
% to memory mapping.  If however, we do require random access to file
% read/write, memory mapping is MUCH faster than fread/fwrite.
if nargin<9||isempty(blocksize)||isequal(blocksize,datasize) % Not stored in blocks
    try
        [elementsize, mmdatatype]=mdatatypeprops(datatype);
        fileobj=memmapfile(filename,'Offset',data_offset,...
            'Format',{mmdatatype [bands_ip*(~isequal(complextype,0)+1) datasize(:).'] 'val'});
        fileobj.data.val(1,1); % Will give out-of-error memory if memmapfile not possible for this data
        chipper_function = @chipper_mm;
        close_function=@() clear('fileobj');
        return;
    end
end
bands_ip=bands_ip*(~isequal(complextype,0)+1); % Complex datatypes have bands for I and Q
% If memmapfile is not possible
clear fileobj; % Release handle on file
fileobj=fopen(filename,'r',endian);
chipper_function = @chipper_fread;
close_function=@() fclose(fileobj);

%% Nested functions
    function data_out=chipper_mm(varargin)
        newargs=reorient_chipper_args(symmetry, datasize, varargin{:});
        data_out=read_bip_mm(fileobj, endian, newargs{:});
        data_out=data2complex(data_out);
        data_out=reorient_chipper_data(symmetry, data_out);
    end

    function data_out=chipper_fread(varargin)
        newargs=reorient_chipper_args(symmetry, datasize, varargin{:});
        if ~exist('blocksize','var')||isempty(blocksize)||isequal(datasize,blocksize) % Band interleaved by pixel
            data_out=read_bip(fileobj, datasize, data_offset, datatype,...
                bands_ip, true, newargs{:});
        else % Band interleaved by block
            data_out=read_bib(fileobj, datasize, data_offset, datatype,...
                bands_ip, blocksize, newargs{:});
        end
        data_out=data2complex(data_out);
        data_out=reorient_chipper_data(symmetry, data_out);
    end

    % Takes data that has been read in in multiple bands and converts it to
    % complex
    function data_out=data2complex(data_in)
        if isequal(complextype,0) % Not complex
            data_out=data_in;
        elseif isequal(complextype,1) % Standard I/Q
            data_out=complex(data_in(:,:,1:2:end),data_in(:,:,2:2:end));
        elseif isa(complextype,'function_handle') % Custom complex types
            data_out=complextype(data_in);
        end
    end


%% Alternate method
% The following commented code would theoretically be a much easier way to
% do build a chipper.  It would leverage code already developed and
% maintained (and theoretically optimized) by MATLAB.  However, MATLAB's
% multibandread function is MUCH, MUCH slower on large datesets (at least
% as of 2010). Until they speed it up, we will use our own code.
%
% chipper_function = @test_chipper;
% close_function = @() 1; % No need to close chipper, so do nothing.
% 
%     function out=test_chipper(varargin)
%         [dim1range,dim2range,subsample]=check_chipper_args(datasize, varargin{:});
%         out=multibandread(filename,...
%             [datasize([2 1]) complexbool+1], datatype, data_offset, 'bip', endian,...
%             {'Row','Range', [dim2range(1), subsample(2), dim2range(2)]},...
%             {'Column','Range', [dim1range(1), subsample(1), dim1range(2) ]});
%         out=permute(out,[2 1 3]);
%         if complexbool
%             out=complex(out(:,:,1:2:end),out(:,:,2:2:end));
%         end
%     end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////