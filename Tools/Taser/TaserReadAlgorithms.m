function Algorithms = ReadAlgorithms()
%ReadAlgorithms Reads Algorithms text file 
%
% ////////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED         ///
% ////////////////////////////////////////////
%
% INPUTS:
%   None  
%
% OUTPUTS:
%  Algorithms - Structure of Algorithms
%
%VERSION:
%   1.0 
%     - Tim Cox 20110613
%     - initial version.
%   1.1
%     - Tim Cox 20130314
%     - Added Override file for compiled users

Algorithms = [];

%determine if user has an Algorithms file in the Temp dir.  This will be
%used to override the standard defaults for a compiled application
OverrideName = [tempdir 'Algorithms.txt'];
if ~isempty(dir(OverrideName))
    filename = OverrideName;
else
    %construct file name
    [pathstr, name, ext] = fileparts(which('TaserReadAlgorithms')); 
    filename = sprintf('%s/Algorithms.txt',pathstr);
end

%determine if file exists
fid = fopen(filename,'r');

%return if fopen was unsucessful
if (fid < 3)
    return;
end

%parse file

%First line is labels
tline = fgetl(fid);

tline = fgetl(fid);
while ischar(tline)
    strings = regexp(tline, ',', 'split');
    temp.TextName = strtrim(strings{1});
    temp.CallName = strtrim(strings{2});
    temp.DataType = strtrim(strings{3});
    temp.Selected = strtrim(strings{4});
    temp.Polarimetric = strtrim(strings{5});
    Algorithms = horzcat(Algorithms,temp);
    tline = fgetl(fid);
end

fclose(fid);

% ////////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED         ///
% ////////////////////////////////////////////