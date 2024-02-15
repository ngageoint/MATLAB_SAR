function test_noise( filename )
%VALIDATE_SICD_NOISE Compare noise description in metadata to actual
%measured noise in pixel data.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('filename','var')
    % Select complex image
    % Recall last interactively selected path used
    if ispref('matlab_sar_toolbox','last_used_directory')
        pathname = getpref('matlab_sar_toolbox', 'last_used_directory');
        if ~ischar(pathname) || ~exist(pathname, 'dir')
            pathname = pwd;
        end
    else
        pathname = pwd;
    end
    [filename, pathname] = uigetfile(sar_file_extensions('complex'), ...
        'Select Input File', pathname, 'MultiSelect', 'off');
    if(filename)
        filename=[pathname filename];
    else % filename=0.  Cancel was pressed, instead of a file being chosen.
        return;
    end
    setpref('matlab_sar_toolbox', 'last_used_directory', pathname);
end

ro = open_reader(filename);
if ~iscell(ro), ro = {ro}; end;
sicd_meta = ro{1}.get_meta();
if ~isfield(sicd_meta, 'Radiometric') || ...
        ~isfield(sicd_meta.Radiometric, 'NoiseLevel') || ...
        ~isfield(sicd_meta.Radiometric.NoiseLevel, 'NoisePoly') || ...
        ~isfield(sicd_meta.Radiometric.NoiseLevel, 'NoiseLevelType') || ...
        ~strcmpi(sicd_meta.Radiometric.NoiseLevel.NoiseLevelType, 'ABSOLUTE')
    file_cleanup();
    error('TEST_NOISE:NO_NOISE_METADATA', 'No noise polynomial found for given data.')
end
mitm_viewer(filename, 'mode', 'aoi', 'selectCallback', @compare_noise);
set(gcf,'DeleteFcn',@file_cleanup); % Close open files on closing
        
    function compare_noise(aoi_info, frame_number)
        azlimits=[aoi_info(1), aoi_info(1)+aoi_info(3)-1];
        rnglimits=[aoi_info(2), aoi_info(2)+aoi_info(4)-1];

        % Multiple frame number returned from mitm_viewer for polarimetric datasets.
        for i = 1:numel(frame_number)
            sicd_meta = ro{frame_number(i)}.get_meta();
            pred_noise = sicd_polyval2d(...
                sicd_meta.Radiometric.NoiseLevel.NoisePoly, ...
                mean(azlimits), mean(rnglimits), sicd_meta);
            disp(['Position: ', num2str(mean(azlimits)), ', ' ...
                num2str(mean(rnglimits)), ' (az, rg)']);
            if numel(frame_number)>1
                disp(['Pol: ', sicd_meta.ImageFormation.TxRcvPolarizationProc]);
            end
            disp(['Predicted noise (dB): ' num2str(pred_noise, 4)]);
            data = ro{frame_number(i)}.read_chip(azlimits, rnglimits);
            sf = sicd_polyval2d(...
                sicd_meta.Radiometric.SigmaZeroSFPoly, ...
                mean(azlimits), mean(rnglimits), sicd_meta);
            meas_noise = 10*log10(sf*sum(double(data(:)).*conj(double(data(:))))/numel(data));
            disp(['Measured noise (dB): ' num2str(meas_noise, 4)]);
        end
    end

    function file_cleanup(obj, eventdata)
        % Clean up intermediate files
        for i=1:length(ro)
            ro{i}.close();
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////