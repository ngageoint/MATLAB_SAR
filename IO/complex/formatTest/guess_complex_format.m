function format_string = guess_complex_format(filename)
% GUESS_COMPLEX_FORMAT Given a file of complex data, this functions guesses
% the format and returns the format type in a string
%
% Assumes a set of test functions (all of the format is*.m) are in the same
% directory.  If format is unrecognized, an empty string is returned.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Get list of all file format test functions
path=fileparts(which(mfilename('fullpath'))); % Path for this mfile
filelist=dir([path filesep 'is*.m']); % All test functions are named with this convention
test_formats = regexprep(sort({filelist.name}),'\.m$',''); % Remove extensions from filenames

% We move the SICD test to the front since 
% 1) It is the foundation of our framework and the format we expect to see
% the most.
% 2) It avoids the regular NITF test triggering.
% 3) Plus the test for non-SICD formats is very quick and light-weight, so
% it does not hurt us when reading non-SICD formats.
sicd_idx = find(strcmp('issicd',test_formats));
test_formats = test_formats([sicd_idx 1:(sicd_idx-1) (sicd_idx+1):end]);

% Run each test until we get positive result
format_string=[];
for i=1:numel(test_formats)
    try
        if(feval(test_formats{i}, filename))
            format_string=test_formats{i}(3:end);
            return
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////