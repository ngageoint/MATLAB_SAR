function demo_core(functionfordemo, varargin)
%DEMO_CORE Call processing algorithm using MATLAB GUIs.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Check for input params
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('MultiSelect','off');
p.FunctionName = 'demo_core';
p.parse(varargin{:});

% Select complex image
% Recall last interactively selected path used
if ispref('matlab_sar_toolbox','last_used_directory')
    pathname = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathname)||~exist(pathname,'dir')
        pathname = pwd;
    end
else
    pathname = pwd;
end
[filename,pathname]=uigetfile(sar_file_extensions('complex'),...
    'Select Input File',pathname,'MultiSelect',p.Results.MultiSelect);
if(iscell(filename)) % Multiple files requested
    for i=1:length(filename), fullfilenames{i}=[pathname filename{i}]; end;
elseif(filename)
    fullfilenames=[pathname filename];
else % filename=0.  Cancel was pressed, instead of a file being chosen.
    return;
end
setpref('matlab_sar_toolbox','last_used_directory',pathname);

mitm_viewer(fullfilenames,'mode','aoi','selectCallback',@calculate_and_display_output);

    function calculate_and_display_output(aoi_info, frame_number)
        azlimits=[aoi_info(1) aoi_info(1)+aoi_info(3)-1];
        rnglimits=[aoi_info(2) aoi_info(2)+aoi_info(4)-1];
        
        % Compute processing
        tmp_name=tempname;
        feval(functionfordemo,fullfilenames,tmp_name,...
            'azlimits',azlimits,'rnglimits',rnglimits,...
            'framenumber',frame_number,varargin{:});
        
        % Get filenames for each frame
        outputfilestructs=dir([tmp_name '*']);
        if(isempty(outputfilestructs)) % Figure was closed without choosing an AOI
            return;
        end
        [resortedname, sortindices]=sort({outputfilestructs.name});
        for j=1:length(outputfilestructs)
            tmp_dir=fileparts(tmp_name);
            fulloutputfilenames{j}=[tmp_dir filesep outputfilestructs(sortindices(j)).name];
        end
        
        % View output
        mitm_viewer(fulloutputfilenames);
        set(gcf,'DeleteFcn',@file_cleanup); % Delete temp files on closing
        
        function file_cleanup(obj, eventdata)
            delete(get(obj,'Children')); % Calls destructor for hg_mitm_viewer, which releases files
            % Clean up intermediate files
            for k=1:length(outputfilestructs)
                delete(fulloutputfilenames{k});
            end
        end
        
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////