function [ element_type, element_length ] = matlabtype2sio( matlabclassstring, iscomplex )
%MATLABTYPE2SIO Convert a MATLAB class string into the element_type and
%element_length fields required for SIO files
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

switch matlabclassstring
    case 'char'
        element_type=1;
        element_length=1;
    case 'single'
        element_type=3;
        element_length=4;
    case 'double'
        element_type=3;
        element_length=8;
    otherwise
        if (strncmpi(matlabclassstring,'uint',4))
            element_type=1;
            element_length=str2double(matlabclassstring(5:end))/8;
        elseif (strncmpi(matlabclassstring,'int',3))
            element_type=2;
            element_length=str2double(matlabclassstring(4:end))/8;
        elseif (strncmpi(matlabclassstring,'float',5))
            element_type=3;
            element_length=str2double(matlabclassstring(6:end))/8;
        else
            error('MATLABTYPE2SIO:UNRECOGNIZED_DATATYPE', 'Reader only recognizes unsigned and signed integers and floats.');
        end
end
if iscomplex % Vector types not supported
    element_type=element_type+10;
    element_length=element_length*2;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////