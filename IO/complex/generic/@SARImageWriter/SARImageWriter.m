% SARImageWriter is an abstract superclass describing the family of all
% writers of SAR raster data.  It was desgined particularly with complex
% SAR data in mind, but subclasses could be built to write detected SAR
% data as well.
% 
% Written by: Wade Schwartzkopf, NGA/ID
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

classdef SARImageWriter < handle
    % Public class data.  Anyone has access to these in a read-only fashion.
    properties(SetAccess=protected, GetAccess=public)
        sicdmeta; % Collection and image formation metadata
        filename;
    end
    
    methods
        % Constructor.  Since this superclass constructor requires
        % arguments, constructors for derived classes will typically call
        % this explicitly (MATLAB's default is to call the supercalss
        % constructor without argument.) to assure that the internal
        % metatdata is correctly set on creation.
        function obj = SARImageWriter(filename, sicdmeta)
            % Fundamental properties
            obj.filename = filename;
            obj.sicdmeta = sicdmeta;
        end
    end
    
    % All subclasses need a way to write pixel data
    methods (Abstract)
        status = write_chip(obj, data, start_indices);
    end
    
    methods (Static)
        % Utility function that mimics the behavior of MATLAB's FSEEK
        % function, but allows for seeking past the end of a file (just
        % appends zeros to the end of the file).
        status = fseek(fid, offset, origin);
    end
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////