function csibatch(varargin)
%CSIBATCH Apply CSI to a set of complex images
%    csibatch('PropertyName',PropertyValue,...)
%
% Displays subaperture type information as color on full resolution data.
%
%       Property name     Description
%       dim               dimension over which to split subaperture
%                            (default = 1)
%       fill              fill factor (default = 1)
%       platformdir       platform direction, 'right' (default) or 'left'.
%                            Assumption is that 2nd dimension is increasing
%                            range.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Select complex image
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end
[filename,pathname]=uigetfile(fullfile(pathstr,'*.*'),'MultiSelect','on');
if iscell(filename)
    number_of_files=length(filename);
else
    if(~filename) % Cancel was chosen
        return;
    else
        number_of_files=1;
        temp=filename; filename=cell(1);
        filename{1}=temp;
    end
end
setpref('matlab_sar_toolbox','last_used_directory',pathname); %store path

h = waitbar(0);
for i=1:number_of_files
    waitbar((i-1)/number_of_files, h, ...
        ['Processing file ' num2str(i) ' of ' num2str(number_of_files) '.']);
    csifile([pathname filename{i}],[pathname filename{i} '.csi'],varargin{:},...
        'azlimits','full','rnglimits','full','rset',true);
end
waitbar(1,h); close(h);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////