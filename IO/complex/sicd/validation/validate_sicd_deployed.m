function validate_sicd_deployed( sicd_filename, report_filename, do_interactive )
%VALIDATE_SICD Generates a validation report for a given SICD
%
% This is a wrapper function that allows you to run validate_sicd.m
% compiled from a command line.
%
% It is recommended that when examining a new source of SICD data, all
% modes from that source be reviewed as well as left- and right-looking
% datasets, if available.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin<3
    do_interactive = true;
else
    do_interactive = logical(eval(do_interactive)); % Convert from string to logical
end

validation_report = validate_sicd( sicd_filename, do_interactive ); % Run validation

% Convert output report in cell array to comma-separate-value (CSV) file
fid = fopen(report_filename, 'w');
for i = 1:size(validation_report,1)
    fprintf(fid,[validation_report{i,1} ', "' validation_report{i,2} '", "' validation_report{i,3} '"\n']);
end
fclose(fid);
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////