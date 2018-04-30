function [ writerobj ] = open_sio_writer( filename, varargin )
%OPEN_SIO_WRITER Initiates a writer for SIO files
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('datatype','single',@(x) ischar(x)); % String from MATLAB class types
% Should datatype default to single, the most common type, or should it
% follow the MATLAB class of the data?
p.addParamValue('complex',[],@(x) isequal(x,1)||isequal(x,0)); % If not specified, use data to determine
p.FunctionName = 'open_sio_writer';
p.parse(varargin{:});

% Variables that should be seen by all subfunctions
% None are set until first write
iscomplex=[];
freadtype=[];
foutid=[]; 
firstfile=[];
dimorder=[];
imsize=[];

writerobj.write_row=@writerow;
writerobj.write_column=@writecolumn;
writerobj.close=@writerclose;

%% Subfunctions
    function writerow(data)
        if isempty(dimorder)
            dimorder=1;
            firstfile=filename;
            setupfile(data);
            imsize=size(data);
            imsize=imsize([2 1]);
        elseif dimorder==1
            if (imsize(2)~=size(data,dimorder))
                error('SIO_WRITE:INVALID_WRITE_SIZE','Row width changed across multiple writes.');
            end
            imsize(1)=imsize(1)+size(data,3-dimorder);
        elseif dimorder==2
            error('SIO_WRITER:INVALID_WRITE_TYPE','Row write after column write not allowed.');
        end
        if iscomplex
            linetosave=zeros(2*numel(data),1,'single');
            linetosave(1:2:end)=real(data);
            linetosave(2:2:end)=imag(data);
        else
            linetosave=data;
        end
        fwrite(foutid,linetosave,freadtype); % Extract angle and write out to file
    end

    function writecolumn(data)
        if isempty(dimorder)
            dimorder=2;
            firstfile=tempname;
            setupfile(data);
            imsize=size(data);
        elseif dimorder==1
            error('SIO_WRITER:INVALID_WRITE','Column write after row write not allowed.');
        elseif dimorder==2
            if (imsize(2)~=size(data,dimorder))
                error('SIO_WRITE:INVALID_WRITE_SIZE','Column height changed across multiple writes.');
            end
            imsize(1)=imsize(1)+size(data,3-dimorder);
        end
        if iscomplex
            linetosave=zeros([2 1].*size(data),'single');
            linetosave(1:2:end)=real(data.');
            linetosave(2:2:end)=imag(data.');
        else
            linetosave=data.';
        end
        fwrite(foutid,linetosave,freadtype); % Extract angle and write out to file
    end

    function writerclose()
        fseek(foutid,4,'bof');
        fwrite(foutid,imsize,'uint32');
        fclose(foutid);
        if(dimorder==2)
            transpose_sio(firstfile,filename);
            delete(firstfile);
        end
    end

    % Setup output file
    function setupfile(data)
        if ~isempty(p.Results.complex)
            iscomplex=p.Results.complex;
        else
            iscomplex=~isreal(data);
        end
        [element_type,element_length]=matlabtype2sio(p.Results.datatype, iscomplex);
        freadtype=p.Results.datatype; % Class type strings for MATLAB classes vary slightly from fread types
        if strcmp(p.Results.datatype,'single')
            freadtype='float32';
        else
            freadtype=p.Results.datatype;
        end
        
        foutid = fopen(firstfile,'w','b');
        fwrite(foutid,hex2dec('FF017FFE'),'uint32'); % SIO magic number
        fwrite(foutid,[0 0],'uint32'); % Write size only at close after we know the size
        fwrite(foutid,element_type,'uint32'); % Element type
        fwrite(foutid,element_length,'uint32'); % Element size (bytes)
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////