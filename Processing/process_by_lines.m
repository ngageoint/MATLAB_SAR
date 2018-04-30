function process_by_lines( compleximagefile, outfile, func_hand, varargin )
%PROCESS_BY_LINES Generic core routine for processing large files in lines
%
% Determines how many lines of the file can be processed in memory at once
% and processes in chunks as large as will fit into memory.
%
%    process_by_lines(compleximagefile, outfile, processing_function_handle, 'PropertyName',PropertyValue,...)
%
%       COMPLEXIMAGEFILE can be either a string that is a filename, or a
%          cell array of open_reader reader objects.
%       OUTFILE is always a string that is a base filename for a set of
%          files (if FUNC_HAND produces a cell array with multiple outputs,
%          multiple output files are produced).
%       FUNC_HAND is a function handle that is applied to every chunk of
%          data.
%
%       Property name     Description
%       azlimits          min and max samples in azimuth over which to
%                            compute (default = query user with GUI).  Use
%                            'full' to specify entire range.  If nothing is
%                            selected, a MATLAB GUI will be used to select
%                            an area.
%       rnglimits         min and max lines in range over which to compute
%                            (default = query user with GUI).  Use 'full'
%                            to specify entire range.  If nothing is
%                            selected, a MATLAB GUI will be used to select
%                            an area.
%       framenumber       frame to process if a multi-image file
%       block_overlap     amount of overlap necessary over adjacent blocks
%                            (default = 0, no overlap)
%       dim               dimension along which processing will be done
%                            (default = 1)
%       function_name     string to display in waitbar while processing
%                            (default = the MATLAB func2str output)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse and validate function arguments
if (nargin<2)||~ischar(outfile)
    error('Output filename required.');
elseif (nargin<1)||(~ischar(compleximagefile)&&~iscell(compleximagefile))
    error('Input filename required.');
end
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('azlimits','default');
p.addParamValue('rnglimits','default');
p.addParamValue('framenumber',1);
p.addParamValue('block_overlap',0,@(x) isvector(x)&&(length(x)<=2));
p.addParamValue('dim',1,@(x) isequal(x,1)||isequal(x,2));
p.addParamValue('function_name',func2str(func_hand),@ischar);
p.addParamValue('outputdatatype','single',@(x) ischar(x)); % String from MATLAB class types
p.FunctionName = 'process_by_lines';
p.parse(varargin{:});
azlimits=p.Results.azlimits;
rnglimits=p.Results.rnglimits;
framenumber=p.Results.framenumber;
if ischar(compleximagefile)
    if(strcmpi(azlimits,'default')&&strcmpi(rnglimits,'default'))
        [aoi_info, framenumber]=mitm_viewer(compleximagefile,'mode','aoi',...
            'closeAfterSelect',true);
        if ~isempty(aoi_info) % User may have closed viewer without selecting AOI
            azlimits=[aoi_info(1) aoi_info(1)+aoi_info(3)-1];
            rnglimits=[aoi_info(2) aoi_info(2)+aoi_info(4)-1];
        end
    end
    fr=open_reader(compleximagefile);
    if iscell(fr) % Close ununsed frames in multiframe data
        ffr{1}=fr{framenumber}; % Save reader for frame of interest
        for i=1:length(fr)
            if i~=framenumber, fr{i}.close(); end;
        end
    else
        ffr{1}=fr;
    end
    clear fr;
    if exist('aoi_info','var')&&isempty(aoi_info)  % AOI GUI started, but AOI not selected
        ffr{1}.close(); clear ffr; return;
    end;
else
    ffr=compleximagefile;
end

metadata=ffr{1}.get_meta(); % Data size is the only metatdata needed.  Assume all images same size.
block_overlap=p.Results.block_overlap(:).'.*[1 1];
dim=p.Results.dim;
otherdim=3-dim; % The non-processing dimension
if(isempty(azlimits)||any(strcmpi(azlimits,{'full','default'})))
    azlimits=[1 double(metadata.ImageData.NumCols)];
end
if(isempty(rnglimits)||any(strcmpi(rnglimits,{'full','default'})))
    rnglimits=[1 double(metadata.ImageData.NumRows)];
end
if ((~verify_limits(azlimits,metadata.ImageData.NumCols))||...
        (~verify_limits(rnglimits,metadata.ImageData.NumRows)))
    error('Invalid subimage range specified.');
end
chipdatasize=[diff(azlimits) diff(rnglimits)]+1;

%% Open input/output files and setup
fw={open_sio_writer(outfile,'datatype',p.Results.outputdatatype)};
if(dim==1)
    fw{1}.writeline=fw{1}.write_row;
else
    fw{1}.writeline=fw{1}.write_column;
end
if(dim==2)
    lines2read=chipdatasize(1);
    firstline=azlimits(1);
    lastline=azlimits(2);
    for i=1:length(ffr)
        read_line{i}=@(lims) ffr{i}.read_chip(lims,rnglimits);
    end
else
    lines2read=chipdatasize(2);
    firstline=rnglimits(1);
    lastline=rnglimits(2);
    for i=1:length(ffr)
        read_line{i}=@(lims) ffr{i}.read_chip(azlimits,lims);
    end
end;
size_to_try=chipdatasize; % First try to read in whole image into memory
size_to_try(otherdim)=size_to_try(otherdim)+(2*block_overlap(otherdim));
lastlineread=firstline-1;

%% Make sure processing block size it at least small enough to handle reading in input data
% If block size is too big for processing, we can adjust size later.
h=waitbar(0,'Calculating amount of memory available...');
error_occurred=1;
while error_occurred
    try
        if is64bit&&(prod(size_to_try)*8>2^25) % Arbitrary size limit for 64-bit systems
            error('Memory chunk too large.');  % Otherwise, systems tries to grab all addressable memory
        end                                    % which can incapacitate a machine
        inputdata=cell(length(ffr),1);
        for i=1:length(inputdata)
            inputdata{i}=complex(zeros(size_to_try(:)','single'));
        end
        current_block=complex(zeros(size_to_try(:)','single'));
        for i=1:3 % Amount of memory needed by typical read_line
            extra{i}=complex(zeros(size_to_try(:)','single'));
        end
        error_occurred=0;
    catch % Then keep trying small amounts until processing is possible
        if size_to_try(otherdim)<(2*block_overlap(otherdim))
            rethrow(lasterror);
        end
        clear extra;
        size_to_try(otherdim)=ceil(size_to_try(otherdim)/2);
    end
end
clear extra;
maxnumlines=size_to_try(otherdim)-(2*block_overlap(otherdim));

%% Progressively process blocks from file with appropriate overlap
% Zeropad beginning and preload buffer so overlap works
new_line_indices=cell(1,2);
new_line_indices{dim}=1:size_to_try(dim);
new_line_indices{otherdim}=(block_overlap(otherdim)+1):2*block_overlap(otherdim);
if block_overlap(otherdim);
    for i=1:length(inputdata)
        inputdata{i}(new_line_indices{:})=read_line{i}(lastlineread+[1 block_overlap(otherdim)]);
    end
end
lastlineread=lastlineread+block_overlap(otherdim);
% Process in blocks with overlap
circshift_indices_to=new_line_indices;
circshift_indices_to{otherdim}=1:2*block_overlap(otherdim);
circshift_indices_from=new_line_indices;
while lastlineread<lastline % For each block
    try % Iterate until we don't get an out of error memory
        numlines2read=min(maxnumlines,lastline-lastlineread);
        new_line_indices{otherdim}=2*block_overlap(otherdim)+(1:numlines2read);
        for i=1:length(inputdata)
            inputdata{i}(new_line_indices{:})=read_line{i}(lastlineread+[1 numlines2read]);
        end
        current_block=func_hand(inputdata{:}); % Process next block
        extra=zeros([size_to_try(:)' 2],'single'); % Test if memory is available for writeline
    catch % Then keep trying to process smaller blocks until processing is possible
        if size_to_try(otherdim)<=(max(1,2*block_overlap(otherdim)))
            rethrow(lasterror);
        end
        size_to_try(otherdim)=ceil(size_to_try(otherdim)/2);
        maxnumlines=size_to_try(otherdim)-(2*block_overlap(otherdim));
        for i=1:length(inputdata)
            oldinputdata=inputdata{i}(circshift_indices_to{:});
            inputdata{i}=complex(zeros(size_to_try(:)','single'));
            inputdata{i}(circshift_indices_to{:})=oldinputdata;
        end
        continue; % Don't write to file if processing failed
    end
    clear extra;
    if(iscell(current_block)) % multiple outputs
        if(length(fw)==1) % Open more writers if necessary
            for k=2:length(current_block)
                fw{k}=open_sio_writer([outfile num2str(k,'%03d')],'datatype',p.Results.outputdatatype);
                if(dim==1)
                    fw{k}.writeline=fw{k}.write_row;
                else
                    fw{k}.writeline=fw{k}.write_column;
                end
            end
        end
        current_blockindices={1:size(current_block{1},1),1:size(current_block{1},2)};
        current_blockindices{otherdim}=block_overlap(otherdim)+(1:numlines2read);
        for i=1:length(current_block)
            fw{i}.writeline(current_block{i}(current_blockindices{:})); % Save to file
        end
    else
        current_blockindices={1:size(current_block,1),1:size(current_block,2)};
        current_blockindices{otherdim}=block_overlap(otherdim)+(1:numlines2read);
        fw{1}.writeline(current_block(current_blockindices{:})); % Save to file
    end
    clear current_block; % read_line needs memory
    circshift_indices_from{otherdim}=numlines2read+(1:2*block_overlap(otherdim));
    for i=1:length(inputdata)
        inputdata{i}(circshift_indices_to{:})=inputdata{i}(circshift_indices_from{:}); % Blocks must overlap for filter to be correct
    end
    lastlineread=lastlineread+numlines2read;
    waitbar(((lastlineread-firstline+1)/lines2read),h,['Calculating ' p.Results.function_name '...']);
end
% We need to zero pad on the end for consistent handling of overlap
if block_overlap(otherdim);
    numlines2read=block_overlap(otherdim);
    zero_lines=2*block_overlap(otherdim)+(1:block_overlap(otherdim));
    for i=1:length(inputdata)
        if(dim==1)
            inputdata{i}(:,zero_lines)=0;
        else
            inputdata{i}(zero_lines,:)=0;
        end
    end
    current_block=func_hand(inputdata{:}); % Process next block
    current_blockindices{otherdim}=block_overlap(otherdim)+(1:numlines2read);
    if(iscell(current_block)) % multiple outputs
        for i=1:length(fw)
            fw{i}.writeline(current_block{i}(current_blockindices{:})); % Save to file
        end
    else
        fw{1}.writeline(current_block(current_blockindices{:})); % Save to file
    end
end

%% Clean up
clear inputdata current_block % Sometimes writer closing requires memory
close(h); % waitbar
for i=1:length(fw)
    fw{i}.close(); % file writer
end
if ischar(compleximagefile)
    ffr{1}.close();
    clear ffr; % Releases memory mapped files
end

end

function range_is_ok = verify_limits(limits, maxsize)
    if(isempty(limits)||strcmpi(limits,'query')||strcmpi(limits,'full'))
        range_is_ok = true;
        return;
    end
    if(isnumeric(limits)&&((isequal(size(limits),[1 2]))||(isequal(size(limits),[1 2]))))
        range_is_ok = (limits(1)>=1)&&(limits(2)<=maxsize);
    else
        range_is_ok = false;
    end
end

function result = is64bit
    compstr=computer;
    result = strcmp(compstr(end-1:end),'64');
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////