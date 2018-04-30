function [ filelist ] = getremaplist( )
%GETREMAPLIST Return a list of all the available remap functions in this
%directory as cell array of strings
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Assumes all remaps are in a subdirectory called 'remap_funcs' under
% the directory containing this script (getremaplist).  Also assumes
% that ONLY remaps functions are in this directory.
path=fileparts(mfilename('fullpath'));
fullfilelist=dir(fullfile(path, 'remap_funcs', '*.m'));

% Strip off the '.m' extension to get the callable
% remap routine name.
filelist=regexprep({fullfilelist.name}, '(.*)\.m', '$1');
end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////