function [ output_meta ] = meta2sicd_xhdr( xhdr )
%META2SICD_XHDR Convert recognized extended header data to sicd metadata
%structure
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

output_meta = struct();
xhdr_fieldnames = fieldnames(xhdr);
for i=1:length(xhdr_fieldnames)
    if exist(strcat('meta2sicd_',xhdr_fieldnames{i},'.m'),'file');
        eval(['output_meta = setstructfields(output_meta, meta2sicd_'...
            xhdr_fieldnames{i} '(xhdr.' xhdr_fieldnames{i} '));']);
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////