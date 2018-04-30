function [ readerobj ] = multisegment_reader(filename, datasize,...
    datatype, complexbool, data_offset, endian, symmetry, bands)
%MULTISEGMENT_READER Creates a reader object for data spread over multiple blocks.
%
% DATASIZE and DATA_OFFSET input variables can both be arrays specifying
% multiple values.  Length of DATA_OFFSET and size(DATASIZE,1) must be
% equal.
%
% Assumes that image is broken up by second dimension (by rows if row-major
% order) and that block info in DATASIZE and DATA_OFFSET are list in order
% for lowest row index to highest.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if length(data_offset)~=size(datasize,1)
    error('MULTISEGMENT_READER:INVALID_INPUT_PARAMETERS',...
        'DATASIZE and DATA_OFFSET must have matching dimensions.');
end

chipper_function=cell(length(data_offset),1);
close_function=cell(length(data_offset),1);
for i=1:length(data_offset)
    [chipper_function{i}, close_function{i}]=generic_chipper( filename,...
        datasize(i,:), datatype, complexbool, data_offset(i), endian,...
        symmetry, bands);
end
rowends=cumsum(datasize(:,2));
rowstarts=[0; rowends(1:end-1)]+1;
fulldatasize=[datasize(1,1) rowends(end)];

readerobj = chipfun2readerobj(@combine_chippers, [datasize(1,1) rowends(end)]);
readerobj.close=@combine_close;


    function output = combine_chippers(varargin)
        [ dim1range, dim2range, subsample ] = check_chipper_args(fulldatasize, varargin{:});
        output=[];
        blocks2chip=[];
        for j=1:length(rowstarts)
            if (dim2range(1)<=rowends(j))&&(dim2range(2)>=rowstarts(j))
                blocks2chip(end+1)=j;
            end
        end
%         showWaitbar=(length(blocks2chip)>1)&&any(subsample>1);
%         if showWaitbar
%             h = waitbar(0);
%         end
        for j=1:length(blocks2chip)
%             if showWaitbar
%                 waitbar((j-1)/length(blocks2chip),h,...
%                     ['Processing block ' num2str(j) ' of ' num2str(length(blocks2chip)) '.']);
%             end
            blockdim2range(1)=max(dim2range(1),rowstarts(blocks2chip(j)));
            blockdim2range(2)=min(dim2range(2),rowends(blocks2chip(j)));
            blockdim2range=blockdim2range-rowstarts(blocks2chip(j))+1;
            output=[output chipper_function{blocks2chip(j)}(dim1range, blockdim2range, subsample)];
        end
%         if showWaitbar, waitbar(1,h); close(h); end;
    end

    function combine_close()
        for j=1:length(close_function)
            close_function{j}();
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////