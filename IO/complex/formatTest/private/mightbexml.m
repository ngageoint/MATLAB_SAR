function [ booleanout ] = mightbexml( filename )
%MIGHTBEXML Test whether a file might be XML or not.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

fid=fopen(filename);
current_char=' ';
while ~feof(fid)&&isspace(current_char)
  current_char=fread(fid,1,'uint8=>char')';
end
fclose(fid);
booleanout=current_char=='<'; % XML file must start with '<' symbol
if booleanout % The first test is extremely quick and filters out most non-XML files
    try % This test is slower, but verifies XML format
        xmlread(filename); % Invalid XML will cause error
    catch % An error will be caught but will still produce red text in command window.  Annoying...
        booleanout=false;
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////