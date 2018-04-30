function transpose_sio( inputfile, outputfile )
%TRANSPOSE_SIO Low-level transpose of SIO formatted data file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid = fopen(inputfile,'r','b');
ihdr = fread(fid,5,'uint32');
trdatasize=ihdr([3 2]);
element_type=ihdr(4);
element_length=ihdr(5);
foutid = fopen(outputfile,'w','b');
fwrite(foutid,hex2dec('FF017FFE'),'uint32'); % SIO magic number
fwrite(foutid,trdatasize(:),'uint32');
fwrite(foutid,element_type,'uint32'); % Element type
fwrite(foutid,element_length,'uint32'); % Element size (bytes)
lastlineread=0;
[matlabclasstype, is_complex, freadtype]=siotype2matlab(element_type,element_length);
complex_factor=is_complex+1;

% Determine how much memory is available to process in the largest blocks possible
error_occurred=1;
size_to_try=trdatasize;
while error_occurred
    try
        out=complex(zeros(size_to_try(:)',matlabclasstype));
        linetosave=zeros([complex_factor 1].*size_to_try(:)',matlabclasstype);
        extra=complex(zeros(size_to_try(:)',matlabclasstype)); % Intermediate computations MATLAB needs
        extra2=complex(zeros(size_to_try(:)',matlabclasstype)); % Intermediate computations MATLAB needs
        error_occurred=0;
    catch
        size_to_try(1)=ceil(size_to_try(1)/2);
    end
end
clear extra*;
maxnumlines=size_to_try(1);

h=waitbar(0,'Transposing data...');
try
    while lastlineread<trdatasize(1)
        numlines2read=min(maxnumlines,trdatasize(1)-lastlineread);
        fseek(fid,lastlineread*element_length+20,'bof');
        out=fread(fid,[complex_factor*numlines2read trdatasize(2)],...
            [num2str(complex_factor*numlines2read) '*' freadtype '=>' matlabclasstype],...
            element_length*(trdatasize(1)-numlines2read));
        if is_complex
            out=complex(out(1:2:end,:),out(2:2:end,:)).';
            linetosave=zeros(2*numel(out),1,matlabclasstype);
            linetosave(1:2:end)=real(out);
            linetosave(2:2:end)=imag(out);
        else
            linetosave=out.';
        end
        
        fwrite(foutid,linetosave,freadtype); % Write out to file
        lastlineread=lastlineread+numlines2read;
        waitbar((lastlineread/trdatasize(1)),h);
    end
catch
    close(h); % Close waitbar
    fclose(fid);
    fclose(foutid);
    rethrow(lasterror); % This is old MATLAB syntax, but we want it to work back to version 6.5
end
close(h);
fclose(fid);
fclose(foutid);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////