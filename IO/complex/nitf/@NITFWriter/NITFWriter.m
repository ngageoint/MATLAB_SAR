% NITFWriter is a class implementing a specific FlatfileImageWriter for
% NITF files using a NITF "container" with a single image segment.  Files
% written as complex will be written as SICD NITFs-- although for SICD
% files, one should really use SICDWriter, which is more robust and handles
% multi-segment SICDs for larger datasets.
%
% LIMITATION: Currently the writer does not handle multi-segment SICDs
% (required for datasets over 10 Gig).
%
% References:
%       MIL-STD-2500C, Department of Defense Interface Standard
%                      National Imagery Transmission Format version 2.1
%                      01 May 2006
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

classdef NITFWriter < FlatfileImageWriter
    
    properties(SetAccess=protected, GetAccess=protected)
        % For now we'll have only one image segment and one data extension
        % segment.  The various headers and subheaders are of fixed sizes.
        NITF_header_length;
        IS_subheader_size;
        image_data_size;
        % DES_header_length  = 200; % This changed between SICD version 0.5 and 1.0.
        DES_header_length  = 973; % Harded-coded from SICD spec
        DES_data;
        total_file_length;
        LUT; % Optional Look Up Table
    end
    
    methods
        function obj = NITFWriter(filename, sicdmeta, varargin)
            % Parse inputs
            % These inputs are really parsed later in the superclass
            % constructor.  Here just need to 1) extract is_complex, which
            % is required for computing header_skip and 2) make sure that
            % no additional arguments are passed-- like header_skip, which
            % could break things, if passed on to the superclass
            % constructor.
            p = inputParser;
            p.addParamValue('data_type',0); % Really parsed later in the superclass constructor
            p.addParamValue('is_complex',true,@(x) isscalar(x)&&islogical(x));
            p.addParamValue('LUT',[]);
            p.FunctionName = 'NITFWriter';
            p.parse(varargin{:});
            
            % Need to compute header_skip before object is created
            NITF_header_length_temp = 417; % Hard-coded to value from NITF spec
            if p.Results.is_complex
                IS_subheader_size_temp  = 512; % Two bands described
            elseif ~isempty(p.Results.LUT)
                IS_subheader_size_temp  = 1272; % One band with 256 element RGB LUT
            else
                IS_subheader_size_temp  = 499; % Only single band described
            end
            header_skip = NITF_header_length_temp + IS_subheader_size_temp;
            
            % Construct object from superclass
            obj = obj@FlatfileImageWriter(filename, sicdmeta, 'header_skip', header_skip, varargin{:});
            
            % Compute NITF-specific properties
            obj.NITF_header_length   = NITF_header_length_temp; % Must be added to object after constructor is called
            obj.IS_subheader_size    = IS_subheader_size_temp; % Must be added to object after constructor is called
            obj.image_data_size      = (obj.is_complex + 1)* ...
                prod(double(obj.image_size)) * ...
                obj.data_size; % Element size
            if obj.image_data_size>9999999999
                error('NITFWRITER:Image_Size_Too_Big','Images larger than 10 Gig not allowed in NITF.');
            end
            obj.DES_data             = sicdstruct2xml(sicdmeta, 'inc_newline',true, ...
                'inc_padding',true, 'pad_depth',3, 'file_type', 'SICD');
            obj.total_file_length    = double(obj.NITF_header_length) + ...
                double(obj.IS_subheader_size) + ...
                double(obj.image_data_size) + ...
                double(obj.DES_header_length) + ...
                double(length(obj.DES_data));
            obj.LUT = p.Results.LUT;
            
            % Open the file and write the NITF header data.  This includes both the
            % NITF header and the Image Segment subheader.  After calling this the
            % file is in a state that data can be written to it (see write_data_chip
            % in the FlatfileImageWriter class).
            obj.FID = fopen(filename,'w', 'b');
            obj.write_nitf_header();
            obj.write_nitf_imsubhdr();
        end
        
        % Destructor
        function delete(obj)
            % Write DES and close file
            offset = obj.NITF_header_length + ...
                obj.IS_subheader_size + ...
                obj.image_data_size;
            obj.fseek(obj.FID, offset, 'bof'); % Seek to end of image data
            obj.write_nitf_dessubhdr();
            fwrite(obj.FID,obj.DES_data);
            
            fclose(obj.FID);
        end
        
    end
    
    methods(Access=protected)
        write_nitf_header(obj);
        lengthWritten = write_nitf_imsubhdr(obj);
        write_nitf_dessubhdr(obj)
        write_nitf_security_tags(obj);
    end
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////