function pdvdemo(varargin)
%PDVDEMO Read image, generate, and view Phase mapping using MATLAB GUIs.
%    pdvdemo('PropertyName',PropertyValue,...)
%    calculates the phase derivative of a complex image (selected through
%    a MATLAB dialog box) using the properties specified. The AOI over
%    which to compute the PDV is selected through a MATLAB GUI.  This
%    version of PDV does NOT require that the complete data fit into
%    memory.  It processes from any format handled by OPEN_READER and
%    outputs to files in SIO format.
%
%       Property name     Description
%       deltax            pixel shift (default = 0.25)
%       filtersize        size of smoothing filter (default = [5 5])
%       filtertype        type of filter, 'mean' (default) or 'median'
%       dim               dimension over which to calculate phase gradient
%                            (default = 1)
%
% Written by: Wade Schwartzkopf
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

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
    'Select Input File',pathname);
if(~filename) % Cancel was chosen
    return;
end
setpref('matlab_sar_toolbox','last_used_directory',pathname);
mitm_viewer([pathname filename],'mode','aoi','selectCallback',@calculate_and_display_output);

    function calculate_and_display_output(aoi_info, frame_number)
        azlimits=[aoi_info(1) aoi_info(1)+aoi_info(3)-1];
        rnglimits=[aoi_info(2) aoi_info(2)+aoi_info(4)-1];
        
        % Compute processing
        tmp_name=tempname;
        tmp_name2=tempname;
        feval(@pdvfile,[pathname filename],tmp_name,... % PDV
            'azlimits',azlimits,'rnglimits',rnglimits,...
            'framenumber',frame_number,varargin{:});
        feval(@chipfile,[pathname filename],tmp_name2,... % Context imagery to flicker
            'azlimits',azlimits,'rnglimits',rnglimits,...
            'framenumber',frame_number,varargin{:});
        
        % View output
        mitm_viewer({tmp_name tmp_name2},'initialRemap','linearremap');
        set(gcf,'DeleteFcn',@file_cleanup); % Delete temp files on closing
        
        function file_cleanup(obj, eventdata)
            delete(get(obj,'Children')); % Calls destructor for hg_mitm_viewer, which releases files
            % Clean up intermediate files
            delete(tmp_name);
            delete(tmp_name2);
        end
        
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////