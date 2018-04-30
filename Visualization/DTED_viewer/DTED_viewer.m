function varargout = DTED_viewer( filename, varargin )
%DTED_VIEWER Elevation data viewer
%   DTED_VIEWER (FILENAME, 'PropertyName', PropertyValue, ...) takes a
%   file of elevation data in DTED format and displays it in a MATLAB 
%   IMSHOW/IMAGESC style viewer
%
%       Property name        Description
%       figureSize           Size of the figure to be drawn. If FIG_SIZE is
%                            not specified, the largest figure size that will
%                            fit in the current screen with an integer zoom
%                            factor is used.
%
%       undulationfilename   The name/path of the file containing WGS84 to
%                            geoid undulations.  This is used in the data
%                            cursor (see geoid_undulation) to display
%                            the height above ellipsoid rather then height
%                            above mean sea level.  If the file doesn't
%                            exist no HOE measurement is shown.
%
%       heightunits          The units to use in the data cursor display of
%                            heights.  Two values are allowed 'feet' and
%                            'meters'.  The default is meters.
%
%
% Author: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Parse input parameters
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('figureSize', [], @(x) isequal(size(x),[1 2])||isequal(size(x),[2 1]));
p.FunctionName = 'DTED_VIEWER';
p.addParamValue('undulationfilename', 'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE');
p.addParamValue('heightunits', 'meters');
p.parse(varargin{:});

%% Open files
%load last path
if ispref('matlab_sar_toolbox','last_used_directory')
    pathstr = getpref('matlab_sar_toolbox','last_used_directory');
    if ~ischar(pathstr)||~exist(pathstr,'dir')
        pathstr = pwd;
    end
else
    pathstr = pwd;
end

[filename,pathname]=uigetfile(fullfile(pathstr,'*.*'),'MultiSelect','off');
    
[elevations, lats, lons, meta] = read_DTED( [pathname filename] );
elevations = flipud(elevations);

setpref('matlab_sar_toolbox','last_used_directory',pathname); %store path

%% Setup initial figure
fig_hand=figure('Name',filename); % Figure handle

im_hand=imagesc(lons,lats,elevations);
colormap('gray');
set(gca,'YDir','normal');

%% Setup figure toolbar callbacks
set(datacursormode(fig_hand),'UpdateFcn',@elevationdatacursorfun);

    function output_txt = elevationdatacursorfun(obj,eventdata)
        if (strcmp(p.Results.heightunits, 'feet'))
            scale = 1/0.3048;
            lab = 'ft';
        else
            scale = 1.0;
            lab = 'm';
        end
        pos        = get(eventdata,'Position');
        elev_index = get(eventdata, 'DataIndex');
        elev_msl   = elevations(elev_index(2),elev_index(1));
        output_txt={[sprintf('Latitude: %.5f deg', pos(2))], ...
                    [sprintf('Longitude: %.5f deg', pos(1))], ...
                    [sprintf('Elevation (MSL): %.0f %s', elev_msl*scale, lab)]};
        if (exist(p.Results.undulationfilename, 'file'))
            undulation = geoid_undulation(pos(2), pos(1), 'undulationFilename', p.Results.undulationfilename);
            elev_hoe   = elev_msl + undulation;
            output_txt={output_txt{:}, [sprintf('Elevation (HAE): %.0f %s', elev_hoe*scale, lab)]};
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////