% SICDWriter is a class implementing a writer for the SICD file format,
% including large multi-segment SICDs.
%
% References:
%       NGA.STND.0024-1_1.0, Sensor Independent Complex Data (SICD), Volume
%          1, Design & Implementation Description Document, 2011-09-28
%       NGA.STND.0024-2_1.0, Sensor Independent Complex Data (SICD), Volume
%          2, File Format Description Document, 2011-08-01
%       MIL-STD-2500C, Department of Defense Interface Standard
%                      National Imagery Transmission Format version 2.1
%                      01 May 2006
%
% Written by: Wade Schwartzkopf, NGA
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

classdef SICDWriter < SARImageWriter
    
    % Constants
    properties(Constant = true, GetAccess=protected)
        ISSIZEMAX = 9999999998; % Image segment size maximum
        ILOCMAX = 99999; % Largest value we can put in image location field

        IS_SUBHEADER_LENGTH  = 512; % Fixed for two bands image segments
        % DES_HEADER_LENGTH  = 200; % Harded-coded from SICD spec (0.5 and before)
        DES_HEADER_LENGTH  = 973; % Harded-coded from SICD spec (1.0)
    end
    
    % All these things are computed at construction time and then will not
    % change after that.
    properties(SetAccess=immutable, GetAccess=protected)
        FID; % Handle to open file
        NITF_header_length = 417; % Default value is for a single segment
        DES_data; % SICD XML string    

        % Image Segment info
        BytesPerRow;
        NumIS; % Number of image segments
        NumRowsIS; % Number of rows in each segment
        FirstRowIS; % Row index of the first row in each segment (zero-based)
        WriterIS; % Writer objects for each segment
    end
    
    methods
        %% Constructor
        function obj = SICDWriter(filename, sicdmeta, varargin)
            % Superclass constructor (must be explicitly called since it requires arguments)
            obj = obj@SARImageWriter(filename, sicdmeta);
            
            % We do not require much from the SICDMETA input parameter.  It
            % does not even have to meet the SICD spec or have valid SICD
            % fields.  At a minimum though, SICDMETA must include the
            % number of rows and columns (ImageData.NumRows and
            % ImageData.NumCols).
            
            % Parse inputs
            % Optional arguments could, in the future, allow for the
            % passing of instructions for adding program specific security
            % tags, TREs, or other NITF components, but we have not
            % implemented that here.
            % p = inputParser;
            % p.addParamValue('optional_argument_name',0);
            % p.FunctionName = 'SICDWriter';
            % p.parse(varargin{:});

            % Compute image segment parameters
            if isfield(sicdmeta, 'ImageData')&&isfield(sicdmeta.ImageData,'PixelType')
                if strcmp(sicdmeta.ImageData.PixelType, 'RE32F_IM32F')
                    BytesPerPixel = 8;
                elseif strcmp(sicdmeta.ImageData.PixelType, 'RE16I_IM16I')
                    BytesPerPixel = 4;
                elseif strcmp(sicdmeta.ImageData.PixelType, 'AMP8I_PHS8I')
                    BytesPerPixel = 2;
                else
                    error('SICDWRITER:UNRECOGNIZED_PIXELTYPE', ...
                        ['PixelType must be RE32F_IM32F, RE16I_IM16I, or AMP8I_PHS8I, but was ' ...
                        sicdmeta.ImageData.PixelType ' instead.']);
                end
            else % If PixelType not found in sicdmeta, default to RE32F_IM32F
                sicdmeta.ImageData.PixelType = 'RE32F_IM32F';
                BytesPerPixel = 8;
            end
            obj.BytesPerRow = double(sicdmeta.ImageData.NumCols) * BytesPerPixel;
            NumRowsLimit = min( floor( obj.ISSIZEMAX / obj.BytesPerRow ), obj.ILOCMAX );
            obj.NumIS = ceil(double(sicdmeta.ImageData.NumRows)/NumRowsLimit);
            obj.FirstRowIS = 0;
            obj.FirstRowIS(2:obj.NumIS) = (1:(obj.NumIS-1)) * NumRowsLimit;
            obj.NumRowsIS(1:(obj.NumIS-1)) = NumRowsLimit;
            obj.NumRowsIS(obj.NumIS) = double(sicdmeta.ImageData.NumRows) - ...
                ((obj.NumIS - 1) * NumRowsLimit);

            % Compute DES parameters
            obj.NITF_header_length = 401 + (16 * obj.NumIS);
            obj.DES_data = sicdstruct2xml(sicdmeta, 'inc_newline', true, ...
                'inc_padding', true, 'pad_depth', 3, 'file_type', 'SICD');
            
            % Open the file and write the NITF file header data.
            obj.FID = fopen(filename,'w', 'b');
            obj.write_sicd_fileheader();

            % Construct writer objects for each segment and link them to
            % our already open file.
            ismeta.ImageData.PixelType = sicdmeta.ImageData.PixelType;
            ismeta.ImageData.NumCols = sicdmeta.ImageData.NumCols;
            for i = 1:obj.NumIS
                ismeta.ImageData.NumRows = obj.NumRowsIS(i);
                obj.WriterIS{i} = FlatfileImageWriter(filename, ismeta, ...
                    'header_skip', obj.NITF_header_length + ...
                    (obj.IS_SUBHEADER_LENGTH*i) + ...
                    (sum(obj.NumRowsIS(1:(i-1))) * obj.BytesPerRow), ...
                    'is_complex', true, ...
                    'fid', obj.FID);
            end
        end
        
        % All of the work done here is distributing the file writing across
        % the multiple NITF image segments in the SICD.  The actual writing
        % to file is done by calling a set of per-segment writer objects
        % that were setup in the constructor.
        function status = write_chip(obj, data, start_indices)
            % We will define firstrows and lastrows of image segments as
            % one-based, rather than zero-based like obj.FirstRowIS.
            % obj.FirstRowIS was defined to reflect the convention in the
            % SICD documents, whereas firstrows and lastrows are defined in
            % such a way to match MATLAB convention and the start_indices
            % input parameter.
            firstrows = obj.FirstRowIS + 1;
            lastrows = firstrows + obj.NumRowsIS - 1;
            
            % Write data to file one segment at a time
            status = 0;
            for j=1:obj.NumIS
                % Is there anything to write in this segment?
                if (start_indices(2) <= lastrows(j)) && ...
                        ((start_indices(2) + size(data,2) - 1) >= firstrows(j))
                    % Indices of rows in entire image that we will be writing
                    rowrange(1) = max(start_indices(2), firstrows(j));
                    rowrange(2) = min(start_indices(2) + size(data,2) - 1, ...
                        lastrows(j));
                    % Indices of rows in DATA input parameter that we will be writing from
                    datarange = rowrange - start_indices(2) + 1;
                    % Indices of NITF image segment that we will be writing to
                    segmentrange = rowrange - firstrows(j) + 1;
                    status = status  + ...
                        obj.WriterIS{j}.write_chip(...
                            data(:,datarange(1):datarange(2)), ...
                            [start_indices(1) segmentrange(1)]);
                end
            end
        end
        
        %% Destructor
        function delete(obj)
            % Write image subheaders upon closing.  We don't do this on
            % open, since if the image is not written yet, jumping to any
            % image subheader beyond the first will result in gigabytes of
            % file being created, which could cause unecessary delay.
            % Another (perhaps better) option would be to write each header
            % the first time any pixel data is written to a segment.
            pos = obj.NITF_header_length;
            for i = 1:obj.NumIS
                obj.fseek(obj.FID, pos,'bof');
                obj.write_sicd_imsubhdr(i);
                pos = pos + obj.IS_SUBHEADER_LENGTH + (obj.NumRowsIS(i) * obj.BytesPerRow);
            end
            % Write DES
            obj.fseek(obj.FID, pos, 'bof'); % Seek to end of image data
            obj.write_sicd_dessubhdr();
            fwrite(obj.FID,obj.DES_data);
            % Close file
            fclose(obj.FID);
        end
        
    end
    
    methods(Access=protected)
        write_sicd_fileheader(obj);
        write_sicd_imsubhdr(obj, n);
        write_sicd_dessubhdr(obj)
        write_sicd_security_tags(obj);
    end
    
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////