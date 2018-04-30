function CropImage(varargin)
% CROPIMAGE Taser wrapper for chipfile.m function
%
% Written by: Tim Cox, NRL
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('filename','');
p.addParamValue('aoi',[0 0 0 0]);
p.addParamValue('segment',1);
p.parse(varargin{:});

filenames = p.Results.filename;
aoi = p.Results.aoi;
segment = p.Results.segment;

%get output folder
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

outdir = uigetdir(pathstr,'Specify Output Folder');

if iscell(filenames)
    for i=1:length(filenames)
        %construct output filename
        [pathstr, name, ext] = fileparts(filenames{i});
        aoistr = sprintf('_%d_%d_%d_%d_seg_%d',aoi(1),aoi(2),aoi(3),aoi(4),segment);
        outname = [outdir filesep name aoistr '.nitf'];

        %crop image
        chipfile(filenames{i},outname,'azlimits',[aoi(1) aoi(1)+aoi(3)-1],'rnglimits',[aoi(2) aoi(2)+aoi(4)-1],...
                 'framenumber',segment,'file_format','SICD','block_size',2^21);
    end
else
    %construct output filename
    [pathstr, name, ext] = fileparts(filenames);
    aoistr = sprintf('_%d_%d_%d_%d_seg_%d',aoi(1),aoi(2),aoi(3),aoi(4),segment);
    outname = [outdir filesep name aoistr '.nitf'];

    %crop image
    chipfile(filenames,outname,'azlimits',[aoi(1) aoi(1)+aoi(3)-1],'rnglimits',[aoi(2) aoi(2)+aoi(4)-1],...
             'framenumber',segment,'file_format','SICD','block_size',2^21);
end

msgbox('Image Crop Completed');

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////