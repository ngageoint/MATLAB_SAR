function format_string = guess_tiff_sensor(tiff_filename)
% GUESS_TIFF_SENSOR Given a TIFF this function tries to look for other
% sensor-specific files in relative paths to the main tiff file that may
% provide more metadata.
%
% Assumes a set of test functions (all of the format istiff*.m) are in the
% same directory.  If format is unrecognized, an empty string is returned.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fullname=which(mfilename('fullpath'));
path=fileparts(fullname);
filelist=dir([path filesep 'istiff*.m']);

format_string=[];
for i=1:length(filelist)
    str_ind=strfind(filelist(i).name,'.m');
    if ~isempty(str_ind)&&(str_ind(end)>1)
        test_format=filelist(i).name(1:(str_ind(end)-1)); % Strip off .m extension
        if(feval(test_format, tiff_filename))
            format_string=test_format(7:end);
            return
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////