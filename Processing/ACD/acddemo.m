function acddemo(early_filename, late_filename)
%ACDDEMO Read image, generate, and view amplitude change detection using MATLAB GUIs.
%    acddemo(early_filename, late_filename)
%
% Interactive GUI for image and AOI selection for amplitude change detection.
%
% Written by: Tom Krauss, NGA/IDT
% Modified by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Select complex images
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end
if nargin<1
    [filename1,pathname1]=uigetfile(fullfile(pathstr,'*.*'),'Select reference image');
    if(~filename1) % Cancel was chosen
        return;
    end
    early_filename=[pathname1 filename1];
    pathstr = pathname1;
    setpref('matlab_sar_toolbox','last_used_directory',pathstr); %store path
end
if nargin<2
    [filename2,pathname2]=uigetfile(fullfile(pathstr,'*.*'),'Select match image');
    if(~filename2) % Cancel was chosen
        return;
    end
    late_filename=[pathname2 filename2];
end

% Select AOI
mitm_viewer(early_filename,'mode','aoi','selectCallback',@calculate_and_display_output);

% Wrapper for actual ACD calculation
    function calculate_and_display_output(aoi_info, frame_number)
        startPoint_1 = [aoi_info(1) aoi_info(2)];
        chipSize = min(2000,aoi_info(3:4));
        skip = [1 1];
        
        full_acd = acdfile(early_filename, late_filename, startPoint_1, chipSize, skip, 'projectToDEM', false);
        full_acd_uint8=uint8(zeros(size(full_acd)));
        for i=1:size(full_acd,3)
            full_acd_uint8(:,:,i)=densityremap(full_acd(:,:,i));
        end
        figure('Name','Slant plane, shadows down'); imshow(permute(full_acd_uint8,[2 1 3]));
        
%         full_acd = acdfile(early_filename, late_filename, startPoint_1, chipSize, skip, 'projectionType', 'ground_northup', 'projectToDEM', false);
%         full_acd_uint8=uint8(zeros(size(full_acd)));
%         for i=1:size(full_acd,3)
%             full_acd_uint8(:,:,i)=densityremap(full_acd(:,:,i));
%         end
%         figure('Name','North up'); imshow(permute(full_acd_uint8,[2 1 3]));
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////