% SIOWriter is a class implementing a specific FlatfileImageWriter for SICD
% files using an SIO "carrier".  The SICD metadata is written in a user
% data field called "SICDMETA" in the standard SIO file format.  The caller
% can also turn off the user data, since not all SIO readers handle this.
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
classdef SIOWriter < FlatfileImageWriter
    
    methods
        % Constructor
        function obj = SIOWriter(filename, sicdmeta, varargin)
            % Parse inputs
            % These inputs are really parsed later in the superclass
            % constructor.  Here just need to 1) extract
            % include_sicd_metadata, which is specific to SIOWriter class
            % and 2) make sure that no additional arguments are passed--
            % like header_skip, which could break things, if passed on to
            % the superclass constructor.
            p = inputParser;
            p.addParamValue('data_type',0);
            p.addParamValue('is_complex',0);
            p.addParamValue('include_sicd_metadata',true,@(x) isscalar(x)&&islogical(x));
            p.FunctionName = 'SIOWriter';
            p.parse(varargin{:});
            
            % We'll use a fixed definition SIO header with (possibly) one user-data segment.
            % The header will the look like this:
            %    Core SIO header:
            %       Magic key   (4 byte uint, fixed value of 'FF027FFD')
            %       Rows        (4 byte uint)
            %       Columns     (4 byte uint)
            %       Data type   (4 byte uint)
            %       Data size   (4 byte uint, # bytes per element)
            %    Optional "user data":
            %       Num pairs   (4 byte uint, # pairs of user data, fixed at 1)
            %       Name bytes  (4 byte uint, # bytes in name of user element,
            %                    fixed at 8 for name "SICDMETA")
            %       Name        (8 bytes containing "SICDMETA")
            %       Value bytes (4 byte uint, value is length of XML string)
            %       Value       (XML string holding SICD metadata)
            header_skip = 5*4; % SIO header with no user data is 5 uint32 words
            if p.Results.include_sicd_metadata
                XML_meta_string = sicdstruct2xml(sicdmeta, 'inc_newline', false, ...
                    'inc_padding', false, 'file_type', 'SICD');
                header_skip = header_skip + 4 + 4 + 8 + 4 + length(XML_meta_string); % Add user data length
            end
            
            % Construct base object
            obj = obj@FlatfileImageWriter(filename, sicdmeta, 'header_skip', header_skip, varargin{:});
            
            % Open file and write SIO header
            obj.FID = fopen(filename,'w', 'b'); % We always write big-endian
            % Magic key
            if p.Results.include_sicd_metadata
                magic_key = hex2dec('FF027FFD'); % Indicates big endian, with user-data
            else
                magic_key = hex2dec('FF017FFE'); % Indicates big endian, with no user-data
            end
            fwrite(obj.FID, magic_key, 'uint32');
            % Rows and columns
            fwrite(obj.FID, obj.image_size(2), 'uint32');
            fwrite(obj.FID, obj.image_size(1), 'uint32');
            % Data type and size
            [ element_type, element_length ] = matlabtype2sio( obj.data_type, obj.is_complex );
            fwrite(obj.FID, element_type, 'uint32');
            fwrite(obj.FID, element_length, 'uint32');
            % User data
            if p.Results.include_sicd_metadata
                fwrite(obj.FID, 1, 'uint32');                       % Num pairs of user data
                fwrite(obj.FID, 8, 'uint32');                       % Pair 1 name length
                fwrite(obj.FID, 'SICDMETA', 'char');                % Pair 1 name
                fwrite(obj.FID, length(XML_meta_string), 'uint32'); % Pair 1 value length
                fwrite(obj.FID, XML_meta_string, 'char');           % Pair 1 value
            end
        end
        
        % Destructor
        function delete(obj)
            % If data has not been written to the end of the pre-defined
            % image data size, fill out rest of file space with zeros.
            image_data_size = (obj.is_complex + 1) * ...
                prod(double(obj.image_size)) * ...
                obj.data_size;
            obj.fseek(obj.FID,obj.header_skip+image_data_size,'bof');

            fclose(obj.FID);
        end
    end
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////