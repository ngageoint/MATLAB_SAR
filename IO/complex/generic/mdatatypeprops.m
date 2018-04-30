function [ elementsize, memmapstr ] = mdatatypeprops( datatypestring )
%MDATATYPEPROPS Returns the elementsize in byte of various MATLAB data
%types used in fread.  Also returns a string valid in memmapfile (which
%isn't always the same).
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

memmapstr=datatypestring;
switch datatypestring
    case 'schar'
        elementsize=1;
        memmapstr='int8';
    case 'uchar'
        elementsize=1;
        memmapstr='uint8';
    case { 'int8', 'uint8' }
        elementsize=1;
    case { 'int16', 'uint16' }
        elementsize=2;
    case { 'int32', 'uint32' }
        elementsize=4;
    case 'float32'
        elementsize=4;
        memmapstr='single';
    case { 'int64', 'uint64', 'float64', 'double' }
        elementsize=8;
        if strcmp(datatypestring,'float64')
            memmapstr='double';
        end
    otherwise
        error('READ_COMPLEX:InvalidArg','Invalid data type.');
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////