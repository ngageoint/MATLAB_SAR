function [ readerobj ] = chipfun2readerobj( chipper_function, datasize )
%CHIPFUN2READEROBJ Converts a chipper functions into a reader object
%
% A chipper function is a simplified function of a file reader.  It takes
% only three arguments: [minIndexInDimension1 maxIndexInDimension1],
% [minIndexInDimension1 maxIndexInDimension1], [subsampleInDimension1,
% subsampleInDimension2].  All other information such as the file handle
% and other auxilliary information must be contained within the chipper
% function.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup local variables
lastrowread=0;
lastcolread=0;

%% Define object methods
readerobj.read_row=@readrow;
readerobj.read_column=@readcolumn;
readerobj.read_chip=@readchip;
readerobj.close=@readerclose;

%% Nested functions
    function out = readrow(linenumber,numlines)
        if (nargin<2)
            numlines=1;
        end
        if ((nargin<1)||isempty(linenumber)) % If LINENUMBER not specified, just reads next line
            linenumber=lastrowread+1;
        end
        out = chipper_function([1 datasize(1)],linenumber+(0:numlines-1),[1 1]);
        lastrowread=linenumber+numlines-1;
    end

    function out = readcolumn(linenumber,numlines)
        if (nargin<2)
            numlines=1;
        end
        if ((nargin<1)||isempty(linenumber)) % If LINENUMBER not specified, just reads next line
            linenumber=lastcolread+1;
        end
        out = chipper_function(linenumber+(0:numlines-1),[1 datasize(2)],[1 1]);
        lastcolread=linenumber+numlines-1;
    end

    function out = readchip(varargin)
        if length(varargin)<4||strncmpi(varargin{4},'none',4)||all(varargin{3}==1)||isempty(varargin{3})
            out = chipper_function(varargin{:});
        else % Apply decimation function to each block
            [dim1range,dim2range,subsample]=check_chipper_args(datasize, varargin{:});
            lengthdim2range=length(dim2range(1):subsample(2):dim2range(2));
            sample_data=chipper_function([dim1range(1) dim1range(1)],...
                [dim2range(1) dim2range(1)]);
            out=zeros(length(dim1range(1):subsample(1):dim1range(2)),...
                lengthdim2range,size(sample_data,3),class(sample_data));
            error_occurred=1;
            blocksize=diff(dim2range)+1; % If whole image works, do it
            maxblocksize=500; % Limit so that a reasonable waitbar interval will show status
            if (blocksize>maxblocksize)
                blocksize=ceil((maxblocksize/2)/subsample(2))*subsample(2); % Assure blocks size is multiple of decimation
                h=waitbar(0,'Reading data...');
            end
            lastlineread=dim2range(1)-1;
            lastlineout=0;
            % Try to process the largest blocks that will fit into memory
            while error_occurred
                try
                    while lastlineread<dim2range(2)
                        if exist('h','var'),
                            waitbar((lastlineread-dim2range(1)+1)/(diff(dim2range)+1),h);
                        end
                        endblockline=min(lastlineread+blocksize,dim2range(2));
                        block_lines = chipper_function([dim1range(1) dim1range(2)],...
                            [lastlineread+1 endblockline]);
                        original_class=class(block_lines);
                        if ~isfloat(block_lines), block_lines=single(block_lines); end; % Many decimation types won't work on ints
                        for i=1:size(block_lines,3)
                            new_size=ceil([size(block_lines,1) size(block_lines,2)]./subsample);
                            block_linescol=im2col(block_lines(:,:,i),subsample,'distinct');
                            block_lines2=feval(varargin{4},block_linescol);
                            out(:,lastlineout+(1:new_size(2)),i)=reshape(block_lines2,new_size);
                        end
                        out=cast(out,original_class);
                        lastlineread=endblockline;
                        lastlineout=lastlineout+new_size(2);
                    end
                    error_occurred=0;
                catch % If the blocksize we attempted was too large, try a smaller one
                    if blocksize<=subsample(2)
                        rethrow(lasterror);
                    end
                    blocksize=ceil((blocksize/2)/subsample(2))*subsample(2); % Assure blocks size is multiple of decimation
                    if ~exist('h','var')
                        h=waitbar(0,'Reading data...'); % Don't display waitbar unless we have to read multiple blocks
                    end
                end
            end
            if exist('h','var'), close(h); end;
        end
    end

    function readerclose()
        % Do nothing on close unless this is overridden
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////