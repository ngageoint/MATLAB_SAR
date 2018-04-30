function Prefs = ReadPreferences()
%ReadPreferences Reads Preferences text file 
%
% ////////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED         ///
% ////////////////////////////////////////////
%
% INPUTS:
%   None  
%
% OUTPUTS:
%  Prefs - Structure of Preferences
%
%VERSION:
%   1.0 
%     - Tim Cox 20110525
%     - initial version.

Prefs = [];

%construct file name
[pathstr, name, ext] = fileparts(which('Taser')); 
filename = sprintf('%s/TaserPreferences.txt',pathstr);

%determine if file exists
fid = fopen(filename,'r');

%return if fopen was unsucessful
if (fid < 3)
    return;
end

%parse file

%ImageOverviewFlag
tline = fgetl(fid);
Prefs.ImageOverview = str2double(tline(strfind(tline,'=')+1:end));
%SICD Location
tline = fgetl(fid);
Prefs.SICDLocation = strtrim(tline(strfind(tline,'=')+1:end));
%Max Resolution
tline = fgetl(fid);
Prefs.MaxResolution = str2double(tline(strfind(tline,'=')+1:end));
%Use Overview
tline = fgetl(fid);
Prefs.UseOverview = str2double(tline(strfind(tline,'=')+1:end));
%PHD NumPulses
tline = fgetl(fid);
Prefs.NumPulses = str2double(tline(strfind(tline,'=')+1:end));
%PHD ProcessingPulses
tline = fgetl(fid);
Prefs.ProcessingPulses = str2double(tline(strfind(tline,'=')+1:end));
%PHD dBCheck
tline = fgetl(fid);
Prefs.dBCheck = str2double(tline(strfind(tline,'=')+1:end));
%Complex Default Template Check
tline = fgetl(fid);
Prefs.ComplexDefaultCheck = str2double(tline(strfind(tline,'=')+1:end));
%Complex Template File
tline = fgetl(fid);
Prefs.ComplexTemplate = strtrim(tline(strfind(tline,'=')+1:end));
%PHD Default Template Check
tline = fgetl(fid);
Prefs.PHDDefaultCheck = str2double(tline(strfind(tline,'=')+1:end));
%PHD Template File
tline = fgetl(fid);
Prefs.PHDTemplate = strtrim(tline(strfind(tline,'=')+1:end));

fclose(fid);

% ////////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED         ///
% ////////////////////////////////////////////