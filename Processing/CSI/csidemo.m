function csidemo(varargin)
%CSIDEMO Read image, generate, and view CSI using MATLAB GUIs.
%    csidemo('PropertyName',PropertyValue,...)
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
% Recall last interactively selected path used
pathname = getpref('matlab_sar_toolbox','last_used_directory', pwd);
if ~ischar(pathname)||~exist(pathname,'dir')
    pathname = pwd;
end
[filename,pathname]=uigetfile(sar_file_extensions('complex'),...
    'Select Input File',pathname);
if(~filename) % Cancel was chosen
    return;
end
setpref('matlab_sar_toolbox','last_used_directory',pathname);
fullfilename=fullfile(pathname,filename);

mitm_viewer(fullfilename,'mode','aoi','selectCallback',@calculate_and_display_output);

    function calculate_and_display_output(aoi_info, frame_number)
        azlimits=[aoi_info(1) aoi_info(1)+aoi_info(3)-1];
        rnglimits=[aoi_info(2) aoi_info(2)+aoi_info(4)-1];
        
        % Compute processing
        tmp_name=tempname;
        feval(@csifile,fullfilename,tmp_name,...
            'azlimits',azlimits,'rnglimits',rnglimits,...
            'framenumber',frame_number,varargin{:});
        
        % Get filenames for each frame
        outputfilestructs=dir([tmp_name '*']);
        if(isempty(outputfilestructs)) % Figure was closed without choosing an AOI
            return;
        end
        for i=1:length(outputfilestructs)
            tmp_dir=fileparts(tmp_name);
            fulloutputfilenames{i}=[tmp_dir filesep outputfilestructs(i).name];
        end
        
        % View output
        mitm_viewer([tmp_name '.mbw']);
        set(gcf,'DeleteFcn',@file_cleanup); % Delete temp files on closing
        
        function file_cleanup(obj, eventdata)
            delete(get(obj,'Children')); % Calls destructor for hg_mitm_viewer, which releases files
            % Clean up intermediate files
            for j=1:length(outputfilestructs)
                delete(fulloutputfilenames{j});
            end
        end
        
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////