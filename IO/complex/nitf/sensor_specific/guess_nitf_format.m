function format_string = guess_nitf_format(imagesubhdr)
% GUESS_NITF_FORMAT Given a NITF image subheader structure, this functions
% tries to determine which set of sensor specific adjustments might need to
% be applied
%
% Assumes a set of test functions (all of the format isnitf*.m) are in the
% same directory.  If format is unrecognized, an empty string is returned.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fullname=which(mfilename('fullpath'));
path=fileparts(fullname);
filelist=dir([path filesep 'isnitf*.m']);

format_string=[];
for i=1:length(filelist)
    str_ind=strfind(filelist(i).name,'.m');
    if ~isempty(str_ind)&&(str_ind(end)>1)
        test_format=filelist(i).name(1:(str_ind(end)-1)); % Strip off .m extension
        if(feval(test_format, imagesubhdr))
            format_string=test_format(7:end);
            return
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////