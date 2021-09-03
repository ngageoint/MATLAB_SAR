function [ output_struct ] = sicdxml2struct( dom_node )
%SICDXML2STRUCT Convert SICD XML into a MATLAB struct
%
% Converts SICD XML into a MATLAB structure that is easy to browse and
% reference.  Makes sure all data types are read as the correct MATLAB
% class, arrays are stored as arrays, etc.
%
% Although this function is called SICDXML2STRUCT, and was originally
% written for SICD XML, it works perfectly well on CPHD XML as well, since
% it is very similar in its structure.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

persistent schema_filename_persistent schema_info; % Store schema info to avoid reparsing

%% Setup XML processing
root_node=dom_node.getDocumentElement;
switch char(root_node.getNodeName)
    % Look for appropriate schema file in path
    case 'CPHD'
        schema_filename = which('CPHD_schema_V1.0.1_2018_05_21.xsd');
    case 'CRSD'
        schema_filename = which('CRSD_schema_V1.0.x_NTB_DRAFT_2021_06_30.xsd');
    case 'SICD'
        schema_filename = which('SICD_schema_V1.1.0_2014_09_30.xsd');
    % case 'SIDD' % Might add this later
    otherwise
        error('SICDXML2STRUCT:UNRECOGNIZED_XML_TYPE','Not a SICD or CPHD XML object.');
end
if ~isempty(schema_filename)
    if isempty(schema_info) || ... % Only parse schema if it hasn't been parsed already
        ~strcmp(schema_filename_persistent, schema_filename)
        % The schema XSD contains data type information for each field in the
        % XML metadata.
        schema_info = parse_sicd_schema(schema_filename);
        schema_filename_persistent = schema_filename;
    end
else
    schema_info = struct('types',[],'master',[]);
end
output_struct = recursfun_xml(root_node, schema_info.master); % Begin processing on root node
if root_node.hasAttribute('xmlns') % Extract and save SICD version
    vers_str=char(root_node.getAttribute('xmlns'));
    vers_str=regexp(vers_str,'(0|[1-9]\d*)(\.(0|[1-9]\d*)){0,3}$','match');
    output_struct.Version=vers_str{1};
    switch char(root_node.getNodeName)
        case 'SICD'
            output_struct=sicd_update_meta(output_struct, output_struct.Version);
        case 'CPHD'
            output_struct=cphd_update_meta(output_struct, output_struct.Version);
    end
end

    %% Recursively process through XML
    function current_struct = recursfun_xml( current_node, schema_struct )
        if current_node.hasAttribute('name') % 'name' attribute requires special handling
            current_struct.name=char(current_node.getAttribute('name')); % Save attribute as subfield
            if any(strcmp(current_node.getNodeName,{'Parameter','Desc'})) % Single text node.  Save as name/value pair
                text_node = current_node.getFirstChild; % Parameter and Desc only allow a single text child
                if isempty(text_node)
                    current_struct.value = '';
                else
                    current_struct.value = char(text_node.getData);
                end
                return; % Value already saved and should be no other children
            else % Has subelements ('GeoInfo','RefPt','ReferencePoint').  Process as usual, just with an extra name field.
            end
        end
        childNodes = current_node.getChildNodes;
        numChildNodes = childNodes.getLength;
        if numChildNodes==0 % Temporary to handle misformed nodes
            current_struct=[];
        end
        for j=1:numChildNodes
            current_child=childNodes.item(j-1);
            current_name=char(current_child.getNodeName);
            switch current_child.getNodeType
                case 1 % Element node
                    % Does schema contain information on this field?
                    if isfield(schema_struct,current_name)
                        child_schema_struct = schema_struct.(current_name);
                        % Find base structure or primitive string
                        while isfield(child_schema_struct, 'SCHEMA_type') && ...
                                isfield(schema_info.types, child_schema_struct.SCHEMA_type)
                            child_schema_struct = schema_info.types.(child_schema_struct.SCHEMA_type);
                        end
                    else
                        % We try to be flexible and read all fields, regardless
                        % of whether the field is described in the schema or
                        % not.  This allows extra custom fields to be included
                        % that may not fit the spec.  Also, if we are reading a
                        % SICD of a different version than the schema we are
                        % using, this should allow that to work (at least
                        % partially) as well.
                        child_schema_struct = [];
                    end
                    if strcmp(current_name,'ICP') % Special case: Index ICP by name, rather than number
                        subfield=char(current_child.getAttribute('index'));
                        subfield=subfield(3:end); % Index values are '1:FRFC', '2:FRLC', '3:LRLC', '4:LRFC'
                        current_struct.(current_name).(subfield)=recursfun_xml(current_child, child_schema_struct);
                    elseif current_child.hasAttribute('index') && ... % Ordered elements
                            ~isempty(current_child.getFirstChild) && ... % Would only happen in a misformed SICD
                            (current_child.getFirstChild.getNodeType == 1 || ... % Not a numeric array
                            isnan(str2double(current_child.getFirstChild.getData)))
                        current_struct.(current_name)(str2double(current_child.getAttribute('index')),1)=...
                            recursfun_xml(current_child, child_schema_struct);
                    elseif (current_child.hasAttribute('exponent1') || ... % Array of numeric values
                            current_child.hasAttribute('index')) && ...
                            ~isempty(current_child.getFirstChild) % Would only happen in a misformed SICD
                        if current_child.hasAttribute('exponent1')
                            index1=str2double(current_child.getAttribute('exponent1'))+1;
                        else
                            index1=str2double(current_child.getAttribute('index'));
                            if strcmp(current_node.getNodeName, 'AmpTable')
                                % Inconsistency in SICD.  The only field in SICD that indexes from 0 rather than 1.
                                index1 = index1 + 1;
                            end
                        end
                        if current_child.hasAttribute('exponent2')
                            index2=str2double(current_child.getAttribute('exponent2'))+1;
                        else
                            index2=1;
                        end
                        current_struct(index1,index2)=recursfun_xml(current_child, child_schema_struct);
                    elseif exist('current_struct','var')&&isfield(current_struct,current_name) % Unordered elements with common name
                        if isfield(child_schema_struct,'Identifier')
                            % CPHD removed has almost no "index" attributes
                            % but a number these named identifiers for
                            % lists.
                            current_struct.(current_name) = [current_struct.(current_name) ...
                                recursfun_xml(current_child, child_schema_struct)];
                        else
                            if ~iscell(current_struct.(current_name)) % Convert to cell array if not already done
                                current_struct.(current_name)={current_struct.(current_name)};
                            end
                            current_struct.(current_name){end+1}=recursfun_xml(current_child, child_schema_struct); % Append to end
                        end
                    else % First (and maybe only) element with this tag
                        current_struct.(current_name)=recursfun_xml(current_child, child_schema_struct);
                    end
                case 3 % Text node
                    if numChildNodes==1
                        % Three ways to get the class of a node
                        if ~isempty(schema_struct) && isfield(schema_struct,'SCHEMA_type') % Current way, from schema
                            class_str = schema_struct.SCHEMA_type;
                        elseif current_node.hasAttribute('class')
                            % Old SICDs (<0.5) used to have class info included in nodes
                            class_str = char(current_node.getAttribute('class'));
                        else % We will have to guess at the class
                            class_str = [];
                        end
                        if ~isempty(class_str) % We know the class type
                            in_string=char(current_child.getData);
                            % XSD datatypes used in CPHD/SICD
                            switch class_str
                                case 'xs:string'
                                    current_struct=in_string; % Do nothing.  Already a string
                                case 'xs:double'
                                    current_struct=str2double(in_string);
                                case {'xs:int','xs:integer'}
                                    current_struct=int32(str2double(in_string));
                                case {'xs:positiveInteger', 'xs:nonNegativeInteger'}
                                    current_struct=uint32(str2double(in_string));
                                case 'xs:dateTime'
                                    % Convert to MATLAB serial date number
                                    temp_str=in_string;
                                    temp_str(isletter(in_string))=' ';
                                    current_struct=datenum(temp_str,31);
                                case 'xs:boolean'
                                    current_struct=eval(lower(in_string)); % 'true' or 'false'
                                otherwise
                                    error('SICDXML2STRUCT:UNKNOWN_CLASS','Unrecognized class type.');
                            end
                        else % Just guessing the class type.  Last resort.
                            current_struct=char(current_child.getData);
                            doubval=str2double(current_struct);
                            if(~isnan(doubval)) % If value is numeric, store as such
                                current_struct=doubval;
                            elseif regexp(current_struct,'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(.)?\d*') % dateTime
                                current_struct(isletter(current_struct))=' ';
                                current_struct=datenum(current_struct,31);
                            elseif any(strcmpi(current_struct,{'true','false'})) % boolean
                                current_struct=eval(lower(current_struct));
                            end
                        end
                    else
                        % Combining text content and children under one
                        % element is discouraged for "data-centric" XML.
                        % We certainly should not find it in SICD XML.
                    end
                otherwise
                    error('SICDXML2STRUCT:INVALID_NODE_TYPE','This type of XML node not expected in SICD.');
            end
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////