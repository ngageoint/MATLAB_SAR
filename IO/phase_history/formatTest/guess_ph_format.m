function format_string = guess_ph_format( filename )
% GUESS_PH_FORMAT Given a file of phase history, this functions guesses the
% format and returns the format type in a string
%
% Assumes a set of test functions (all of the format is*.m) are in the same
% directory.  If format is unrecognized, an empty string is returned.
%
% Often file formats for phase history are not as simple to identify as
% their complex data counterparts.  Sometimes filename patterns must be
% used and in some cases the phase history is stored as a set of files in a
% directory (in this case, a single file in that directory should be
% selected).  For this reason, this format guessing function is not quite
% as robust as GUESS_COMPLEX_FORMAT, but it seems to work for the limited
% formats we have.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fullname=which(mfilename('fullpath'));
path=fileparts(fullname);
filelist=dir([path filesep 'is*.m']);

format_string=[];
for i=1:length(filelist)
    str_ind=strfind(filelist(i).name,'.m');
    if ~isempty(str_ind)&&(str_ind(end)>1)
        test_format=filelist(i).name(1:(str_ind(end)-1)); % Strip off .m extension
        if(feval(test_format, filename))
            format_string=test_format(3:end);
            return
        end
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////