function [ matlabclasstype, iscomplex, freadtype ] = siotype2matlab( element_type, element_length )
%SIOTYPE2MATLAB Convert the SIO description of data type in a MATLAB class
%string
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

iscomplex=floor(element_type/10); % Tens place is zero => real, tens place is one => complex
if iscomplex>1  % Tens place is two => vector
    error('OPEN_SIO_READER:VECTOR_TYPE','Vector types for SIO files are not supported by this reader.');
end
datatypenum=mod(element_type,10); % 1: unsigned int, 2: signed int, 3: float
switch datatypenum
    case 1
        datatype = 'uint';
    case 2
        datatype = 'int';
    case 3
        datatype = 'float';
    otherwise
        error('OPEN_SIO_READER:UNRECOGNIZED_DATATYPE', 'Reader only recognizes unsigned and signed integers and floats.');
end
datalength=element_length*8;
if(iscomplex), datalength=datalength/2; end;
freadtype=[datatype num2str(datalength)];

matlabclasstype=freadtype;
switch freadtype
    case 'float32'
        matlabclasstype='single';
    case 'float64'
        matlabclasstype='double';
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////