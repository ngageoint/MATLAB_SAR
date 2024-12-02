function [ schema_struct ] = parse_sicd_schema( filename )
% PARSE_SICD_SCHEMA Parse SICD/CPHD schema XSD into a MATLAB structure
% It is MUCH faster to traverse through a MATLAB structure than DOM nodes,
% so we want to convert the schema info into a MATLAB structure before
% traversing through the SICD/CPHD XML.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Open file
schema_root = xmlread(filename);
schema_root = schema_root.getDocumentElement;
schema_struct.ns = char(schema_root.getAttribute('xmlns'));

% Separate type definitions and master structure definition
for i=1:schema_root.getLength
    child=schema_root.item(i-1);
    if child.getNodeType == child.ELEMENT_NODE
        switch char(child.getNodeName)
            case {'xs:simpleType','xs:complexType'} % Type definitions
                schema_struct.types.(char(child.getAttribute('name'))) = recursfun_schema(child);
            case 'xs:group'
                schema_struct.groups.(char(child.getAttribute('name'))) = recursfun_schema(child);
            case 'xs:element' % Master node (should be onle one)
                schema_struct.master = recursfun_schema(child);
            otherwise
                % Should not be in other types of nodes in top level
                error('SICDXML2STRUCT:UNEXPECTED_NODE_TYPE','Unexpected node type in top level of XSD.');
        end
    end
end

if ~isfield(schema_struct, "groups")
  schema_struct.groups = [];
end
schema_struct.types = recursfun_group(schema_struct.types, schema_struct.groups);

    function output_struct = recursfun_schema(current_node)
        output_struct = struct();
        for j=1:current_node.getLength
            current_child=current_node.item(j-1);
            if current_child.getNodeType == current_child.ELEMENT_NODE
                switch char(current_child.getNodeName)
                    case 'xs:element'
                        name_str = char(current_child.getAttribute('name'));
                        type_str = char(current_child.getAttribute('type'));
                        if ~isempty(type_str)
                            output_struct.(name_str).SCHEMA_type = type_str;
                        else % Element with empty type.  Should have a structure defined within it.
                            output_struct.(name_str) = recursfun_schema(current_child);
                        end
                    case {'xs:restriction','xs:extension'}
                        output_struct = recursfun_schema(current_child); % Adds any attributes
                        output_struct.SCHEMA_type = char(current_child.getAttribute('base'));
                    case {'xs:simpleType','xs:simpleContent',...
                            'xs:complexType','xs:complexContent'}
                        output_struct = recursfun_schema(current_child);
                    case {'xs:group'}
                        ref_str = char(current_child.getAttribute('ref'));
                        output_struct.(ref_str).SCHEMA_group = ref_str;
                    case {'xs:sequence','xs:choice','xs:all'}
                        output_struct = setstructfields(output_struct, recursfun_schema(current_child));
                    case {'xs:attribute'}
                        if isfield(output_struct,'SCHEMA_attributes')
                            output_struct.SCHEMA_attributes = ...
                                [output_struct.SCHEMA_attributes, ...
                                char(current_child.getAttribute('name'))];
                        else
                            output_struct.SCHEMA_attributes = ...
                                {char(current_child.getAttribute('name'))};
                        end
                    case {'xs:minInclusive','xs:maxInclusive',...
                            'xs:minExclusive','xs:maxExclusive',...
                            'xs:enumeration','xs:pattern'}
                        % These fields are expected, but we don't use them
                        % for anything.
                    otherwise
                        error('SICDXML2STRUCT:UNRECOGNIZED_NODE_TYPE','Unrecognized node type in XSD.');
                end
            end
        end
    end
end

function current_struct = recursfun_group(current_struct, groups)
    fields = fieldnames(current_struct);
    for f = 1:numel(fields)
        current_field = current_struct.(fields{f});
        if isfield(current_field, "SCHEMA_group")
            current_struct = assign_group(current_struct, groups);
            current_struct = recursfun_group(current_struct);
            break;
        elseif isfield(current_field, "SCHEMA_type") || isfield(current_field, "SCHEMA_attributes")
            continue;
        elseif ~isempty(fieldnames(current_field))
            current_struct.(fields{f}) = recursfun_group(current_field, groups);
        else
            error('SICDXML2STRUCT:UNEXPECTED_SCHEMA_TYPE','Unexpected schema node type in XSD after initial parsing.')
        end        
    end
end

function new_struct = assign_group(current_struct, groups)
    fields = fieldnames(current_struct);
    for f = 1:numel(fields)
        current_field = current_struct.(fields{f});
        if isfield(current_field, "SCHEMA_group")
            current_group = groups.(current_field.SCHEMA_group);
            group_fields = fieldnames(current_group);
            for gf = 1:numel(group_fields)
              new_struct.(group_fields{gf}) = current_group.(group_fields{gf});
            end
        else
            new_struct.(fields{f}) = current_field;
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////