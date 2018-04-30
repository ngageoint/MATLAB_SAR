function write_sio(filename, data, varargin)
%WRITE_SIO Write SIO formatted complex data
%   Very simple SIO writer.  Only allows a single array held in memory to
%   be written.  Also requires file be big-endian and does not allow for
%   user data.  For more sophisticated SIO writing, use SIOWriter.
%
%   complexdata = write_sio( filename, data )
%
%   INPUTS:
%      FILENAME:      Input file to read.
%      DATA:          Complex array of data to write to SIO file
%      ELEMENTTYPE:   MATLAB data type to write to SIO file (default
%                        'single')
%
%   Author: Wade Schwartzkopf (NGA/IDT)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if(nargin>2)
    datatype=varargin{1};
else
    datatype='single';
end

iscomplex=~isreal(data);
[element_type,element_length]=matlabtype2sio(datatype, iscomplex);
freadtype=datatype; % Class type strings for MATLAB classes vary slightly from fread types
if strcmp(datatype,'single')
    freadtype='float32';
end

%write, big-endian, all non-header data 32-bit float
fid = fopen(filename,'w','b');
fwrite(fid,hex2dec('FF017FFE'),'uint32'); % SIO magic number
fwrite(fid,[size(data,2) size(data,1)],'uint32'); % Write size
fwrite(fid,element_type,'uint32'); % Element type
fwrite(fid,element_length,'uint32'); % Element size (bytes)

if iscomplex
    linetosave=zeros(2*numel(data),1,'single');
    linetosave(1:2:end)=real(data);
    linetosave(2:2:end)=imag(data);
else
    linetosave=data;
end
fwrite(fid,linetosave,freadtype); % Write data to file

fclose(fid);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////