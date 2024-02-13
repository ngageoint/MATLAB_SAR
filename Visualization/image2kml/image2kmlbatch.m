function image2kmlbatch(imagedirname, outputfilename, varargin)
%IMAGE2KMLBATCH Generates KML output decribing all images in a specified directory.
%
% Example call:
%    image2kmlbatch('c:\junk', 'c:\junk\testing2.kml', ...
%                   'name', 'testGoogleEarthLabel', ...
%                   'filter_function', @example_graze_filt)
%
%
% Inputs:
%       imagedirname:   The name of the directory in which to search for images
%       outputfilename: The full path of the resulting KML file
%       PropertyName:   The PropertyName and PropertyValue are key-value
%       PropertyValue:  pairs used throughout the tool kit...
%
% Allowed properties:
%       Property name         Description
%       -----------------------------------------------------------------
%       name            The name to put into the KML file.  This name is
%                       displayed in the "Places" tree viewer by Google
%                       Earth.
%       filter_function The function handle pointing to the SICD metadata
%                       filter function.  This should point to a function
%                       that takes one input (or two if searching for
%                       paired collects), SICD structure(s), and returns a
%                       boolean.  A return of 'true' indicates the image(s)
%                       should be included in the KML generation, false
%                       means "don't include it."
%       In addition, all properties accepted by add_sar_2kml are allowed,
%       except overlay_filename since there will be multiple.
%
% Written by: Tom Krauss, NGA/IDT; Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('filter_function', @(x) true);
p.addParamValue('name', 'Untitled');
p.FunctionName = mfilename;
p.parse(varargin{:});

% Open the directory and load the metadata from each image found.  The
% returned value is a cell array of structures
if nargin(p.Results.filter_function)==1
    image_data = image_search(imagedirname, p.Results.filter_function);
elseif nargin(p.Results.filter_function)==2
    image_data = image_pair_search(imagedirname, p.Results.filter_function);
    pair_filenames(1:2:2*length(image_data))={image_data.filename1};
    pair_filenames(2:2:2*length(image_data))={image_data.filename2};
    pair_meta(1:2:2*length(image_data))={image_data.meta1};
    pair_meta(2:2:2*length(image_data))={image_data.meta2};
    image_data=struct('filename',pair_filenames,'meta',pair_meta);
end
% fprintf('Found %d images\n', length(image_data));

k = kml(p.Results.name); %create an kmltoolbox object
k.filename = outputfilename;
attempt_overlays = isfield(p.Unmatched,'overlay_max_size') && ...
    (p.Unmatched.overlay_max_size>0);
k.zip = attempt_overlays; % Only zip if we are attempting to create overlay images

for i=1:length(image_data)
    if attempt_overlays
        new_args = p.Unmatched;
        new_args.overlay_filename = sprintf('Overlay_%d.png',i); % Will be put in same directory where KML is
        add_sar_2kml(k,fullfile(imagedirname,image_data(i).filename),new_args);
    elseif isfield(image_data(i).meta,'PVP') % Phase history file
        add_sar_2kml(k,fullfile(imagedirname,image_data(i).filename),p.Unmatched);
    else
        add_sar_2kml(k,image_data(i).meta,p.Unmatched);
    end
end

k.save();
if ~isempty(k.includeFiles)
    delete(k.includeFiles{:});
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////