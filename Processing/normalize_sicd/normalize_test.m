function normalize_test( filename )
%NORMALIZE_TEST Function to test complex data normalization
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Select file if not passed in
if ~exist('filename','var')
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
    if ~filename, return; end; % Cancel was pressed
    filename = fullfile(pathname,filename);
    setpref('matlab_sar_toolbox','last_used_directory',pathname);
end

% Select AOI and pass to normalization function
mitm_viewer(filename,'mode','aoi','selectCallback',@normalization_analysis);

    function normalization_analysis(aoi_info, frame_number)
        %% Compute normalized data
        azlimits = [aoi_info(1) aoi_info(1)+aoi_info(3)-1];
        rnglimits = [aoi_info(2) aoi_info(2)+aoi_info(4)-1];
        complex_data = single(read_complex_data(filename, azlimits, rnglimits, ...
            [1 1],'none',frame_number));
        deskew_data1 = deskewfile(filename, 'dim', 1, ...
            'azlimits', azlimits, 'rnglimits', rnglimits, ...
            'framenumber', frame_number );
        normalized_data1 = normalize_complex_file(filename, 'dim', 1, ...
            'azlimits', azlimits, 'rnglimits', rnglimits, ...
            'framenumber', frame_number );
        deskew_data2 = deskewfile(filename, 'dim', 2, ...
            'azlimits', azlimits, 'rnglimits', rnglimits, ...
            'framenumber', frame_number );
        normalized_data2 = normalize_complex_file(filename, 'dim', 2, ...
            'azlimits', azlimits, 'rnglimits', rnglimits, ...
            'framenumber', frame_number );

        %% Azimuth normalization
        % Show FFT before and after
        figure('Name','Azimuth','MenuBar','none','Toolbar','none');
        subplot(2,2,1); imshow(abs(fftshift(fft2(complex_data))).',[]);
        title('Before Normalization');
        subplot(2,2,2); imshow(abs(fftshift(fft2(normalized_data1))).',[]);
        title('After Normalization');
        % Show estimated weighting before and after
        subplot(2,2,3); plot(fftshift(sum(abs(fft(deskew_data1,[],1)),2)));
        title('Before Deweighting');
        subplot(2,2,4); plot(fftshift(sum(abs(fft(normalized_data1,[],1)),2)));
        title('After Deweighting');

        %% Range normalization
        % Show FFT before and after
        figure('Name','Range','MenuBar','none','Toolbar','none');
        subplot(2,2,1); imshow(abs(fftshift(fft2(complex_data))).',[]);
        title('Before Normalization');
        subplot(2,2,2); imshow(abs(fftshift(fft2(normalized_data2))).',[]);
        title('After Normalization');
        % Show estimated weighting before and after
        subplot(2,2,3); plot(fftshift(sum(abs(fft(deskew_data2,[],2)),1)));
        title('Before Deweighting');
        subplot(2,2,4); plot(fftshift(sum(abs(fft(normalized_data2,[],2)),1)));
        title('After Deweighting');
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////