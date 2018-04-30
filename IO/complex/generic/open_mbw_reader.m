function [ readerobj ] = open_mbw_reader(filename)
%OPEN_MBW_READER Intiates a reader object for a multi-band wrapper file
%
% An intermediate format just for handling data that is contained in
% multiple files (color, polarimetric, etc.)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Read mbw fields
fid=fopen(filename,'r');
while ~feof(fid)
    [fieldname,value]=strtok(fgetl(fid),'=');
    if ~isempty(value)
        mbw.(fieldname)=value(2:end);
        doubval=str2double(mbw.(fieldname));
        if(~isnan(doubval)) % If value is numeric, store as such
            mbw.(fieldname)=doubval;
        end
    end
end
fclose(fid);
if(strcmp(mbw.TYPE,'RGB'))
    for i=1:3, filenamelist{i}=mbw.(mbw.TYPE(i)); end
    readerobj=open_stacked_set(filenamelist{:});
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////