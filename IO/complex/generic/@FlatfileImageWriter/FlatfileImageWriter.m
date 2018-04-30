% FlatfileImageWriter is a class describing a family of format-specific
% image writers.  It describes data that can be written to a "flat" file
% (continuous pixels stored in raster order on disk), possibly with some
% header and footer, where data is written with the MATLAB fwrite function.
%
% The FID property of this class assumes an fread/fwrite functionality.  If
% memory-mapped writing, or another library for writing is used, FID (and
% possibly a number of the other properties) would not be appropriate.
% In this case a more generic superclass (SARImageWriter) should be used.
%
% Derived classes should have a constructor that calls superclass
% constructor with appropriate values, opens file, writes header, etc. and
% a destructor that writes any footer info, closes file, etc.
% 
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

classdef FlatfileImageWriter < SARImageWriter
    
    % Protected class data.  Accessed by member classes.
    properties(SetAccess=protected, GetAccess=protected)
        FID;
    end
    
    % Private/protected class data.  No public acess to these.  Their values
    % are set on initialization.
    properties(SetAccess=private, GetAccess=protected)
        % Input parameters that are required to define the structure of the
        % "flat" file to write.
        header_skip;       % Number of bytes to skip before writing data
        data_type;         % "uint8", "float32", etc.
        is_complex;        % Boolean
        
        % Derived parameters computed for ease of use
        image_size;        % Size of image in columns, rows (2-vector), derived from sicdmeta NumCols, NumRows
        data_size;         % Size of image data in bytes
        buffer_type;       % MATLAB data class of DATA_TYPE (so 'float32' becomes 'single'...)
    end
    
    methods
        % Constructor.  Constructors for derived classes will typically
        % call this explicitly to assure that the internal metatdata is
        % correctly set on creation.
        function obj = FlatfileImageWriter(filename, sicdmeta, varargin)
            % Superclass constructor (must be explicitly called since it requires arguments)
            obj = obj@SARImageWriter(filename, sicdmeta);
            
            % Parse constructor parameters
            p = inputParser;
            p.KeepUnmatched = true; % Some subclasses may have extra args which are passed through, which we will ignore
            p.addParamValue('header_skip',0,@(x) isscalar(x));
            p.addParamValue('data_type',get_default_data_type(sicdmeta),@(x) ischar(x));
            p.addParamValue('is_complex',true,@(x) isscalar(x)&&islogical(x));
            p.addParamValue('fid',[],@(x) isscalar(x)); % Points to an already open file.  Although perhaps this case should be a subclass...
            p.FunctionName = 'FlatfileImageWriter';
            p.parse(varargin{:});
            
            % Fundamental properties
            obj.header_skip = p.Results.header_skip;
            obj.data_type = p.Results.data_type;
            obj.is_complex = p.Results.is_complex;
            obj.FID = p.Results.fid;
            
            % Derived properties
            % SICDMETA must have, at a minimum, the number of rows and columns
            obj.image_size  = double([sicdmeta.ImageData.NumCols sicdmeta.ImageData.NumRows]);
            [ obj.data_size, obj.buffer_type ] = mdatatypeprops(obj.data_type);
            
            % Helper function to be used for parameter validation
            function default_data_type = get_default_data_type(sicdmeta)
                default_data_type = 'float32'; % If not passed in as parameter or found in sicdmeta, use this
                if isfield(sicdmeta, 'ImageData')&&isfield(sicdmeta.ImageData,'PixelType')
                    if strcmp(sicdmeta.ImageData.PixelType, 'RE16I_IM16I')
                        default_data_type = 'int16';
                    elseif strcmp(sicdmeta.ImageData.PixelType, 'RE32F_IM32F')
                        default_data_type = 'float32';
                    elseif strcmp(sicdmeta.ImageData.PixelType, 'AMP8I_PHS8I')
                        default_data_type = 'uint8';
                    else
                        % Error: Unrecognized data type
                    end
                end
            end
        end
        
        % Generic data chip writing routine.  Use this to write chips to an
        % already-opened file.  Note that the file is opened with a call to
        % the (abstract) open function.  The file handle is then stored with
        % the object so it doesn't get passed around with this routine.
        status = write_chip(obj, data, start_indices);
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////