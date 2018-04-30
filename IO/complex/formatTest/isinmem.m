function boolout = isinmem(filename)
% Most format test functions refer to a file.  This function test the
% special case where a user passes an array of data in memory, rather than
% through a file.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

boolout = isnumeric(filename) || ... % Just a bare array
    (iscell(filename) && numel(filename)>1 && ... % Cell array
    isnumeric(filename{1}) && isstruct(filename{2}) && ... % with data array and SICD metadata structure
    isfield(filename{2},'ImageData') && ... % Check for minimum required SICD fields
    all(isfield(filename{2}.ImageData,{'NumRows','NumCols'})));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////