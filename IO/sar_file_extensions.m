function [ file_ext_cell_array ] = sar_file_extensions( data_types )
%SAR_FILE_EXTENSIONS Returns a set of valid file extensions for SAR data
%files in a cell array that can be used in MATLAB's UIGETFILE function
%
% Input DATA_TYPES is a string or a cell array of strings with the types of
% SAR data desired.  Options can be 'complex', 'phd' (phase history data),
% or a cell array with both.  Default is 'complex'.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin==0
    data_types = 'complex';
end

file_ext_cell_array = {};

% Complex data file formats
if any(strcmpi('complex',data_types))
    complex_types_cell = {...
        '*.ntf;*.nitf;', 'SICD/NITF (*.ntf,*.nitf)';...
        '*.1__A;*.1__D;', 'ALOS-PALSAR 2 (*.1__A,*.1__D)';...
        '*.h5;', 'COSMO-SkyMed (*.h5)';...
        '*.gff;', 'GFF (*.gff)';...
        '*.tif;*.xml;', 'RadarSat-2 (*.tif;*.xml)';...
        '*.tiff;*.safe;', 'Sentinel-1 (*.tiff;*.safe)';...
        '*.sio;', 'SIO (*.sio)';...
        '*.cos;*.xml;', 'TerraSAR-X COSAR (*.cos,*.xml)';...
        };
    complex_all_extensions = {[complex_types_cell{:,1}], 'SAR Complex Data'};
    file_ext_cell_array = [file_ext_cell_array; complex_all_extensions];
end

% Phase history file formats
if any(strcmpi('phd',data_types))
    phd_types_cell = {...
        '*.cph;*.cphd;', 'CPHD (*.cph,*.cphd)';...
        '*.crsd;', 'CRSD (*.crsd)'...
        };
    phd_all_extensions = {[phd_types_cell{:,1}], 'SAR Phase History'};
    file_ext_cell_array = [file_ext_cell_array; phd_all_extensions];
end

% Add individual types
if any(strcmpi('complex',data_types))
    file_ext_cell_array = [file_ext_cell_array; complex_types_cell];
end
if any(strcmpi('phd',data_types))
    file_ext_cell_array = [file_ext_cell_array; phd_types_cell];    
end

% Add all types for both complex and phase history
if any(strcmpi('complex',data_types))&&any(strcmpi('phd',data_types))
    all_types = [[complex_types_cell{:,1}] [phd_types_cell{:,1}]];
    file_ext_cell_array = [ {all_types, 'SAR Complex Data, Phase History'};...
        file_ext_cell_array]; 
end

% Add "All files" option
file_ext_cell_array = [file_ext_cell_array; {'*.*','All Files (*.*)'}];

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////