function xmlstr = sicdstruct2xml(sicdmeta,varargin)
%SICDSTRUCT2XML Convert SICD or CPHD MATLAB struct into an XML string
%
% Converts SICD MATLAB structure as created by the open_reader framework
% into an XML string that can be written to a SICD file.
%
% Although this function is called SICDSTRUCT2XML, and was originally
% written for SICD XML, it works perfectly well on CPHD XML as well, since
% it is very similar in its structure.
%
% Author: Wade Schwartzkopf NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Parse input parameters
p1 = inputParser;
p1.KeepUnmatched=true;
p1.addParamValue('inc_newline',true);
p1.addParamValue('inc_class_attributes', false, @(x) isscalar(x)&&islogical(x));
p1.addParamValue('file_type','SICD',@(x) any(strcmp(x,{'SICD','CPHD','CRSD'})));
p1.parse(varargin{:});

% Create XML
doc = com.mathworks.xml.XMLUtils.createDocument(p1.Results.file_type);
root_node = doc.getDocumentElement;
switch p1.Results.file_type
    case 'SICD'
        schema_filename = which('SICD_schema_V1.2.1_2018_12_13.xsd');
        % SICD version also written in write_nitf_dessubhdr.m
    case 'CPHD'
        schema_filename = which('CPHD_schema_V1.0.1_2018_05_21.xsd');
    case 'CRSD'
        schema_filename = which('CRSD_schema_V1.0.x_NTB_DRAFT_2021_06_30.xsd');
end
if ~isempty(schema_filename)
    % The schema XSD contains data type information for each field in the
    % XML metadata.
    schema_info = parse_sicd_schema(schema_filename);
else
    schema_info = struct('ns',[],'types',[],'master',[]);
end
if ~isempty(schema_info.ns)
    root_node.setAttribute('xmlns',schema_info.ns);
end
sicdstruct2xml_recurse(root_node, sicdmeta, schema_info.master);
xmlstr = xmlwrite(root_node);
% Remove XML declaration which xmlwrite adds unavoidably:
% <?xml version="1.0" encoding="utf-8"?>
start_index = regexp(xmlstr,'(?<=^<\?xml.+?\?>[\s]*)<');
% Regular expression explanation:
% (?<=    Start of lookbehind
% ^       Start of line (All XML declarations must start first position of first line.)
% <\?xml  Start of XML declaration must be '<?xml'
% .+?     Lazy match for anything with XML declaration
% \?>     End of XML declaration must be '?>'
% [\s]*   Ignore whitepaces immediately after XML declaration
% )       End of lookbehind
% <       We want the beginning of the XML after the XML declaration
if ~isempty(start_index)
    xmlstr = xmlstr(start_index(1):end);
end
% The Saxon XML processor that MATLAB's xmlwrite uses automatically adds
% newlines and tabs.
if ~p1.Results.inc_newline
    xmlstr = regexprep(xmlstr,'[\n][\s]*',''); % Remove newlines and whitespace if requested
end

% Validate XML string against schema
if ~isempty(schema_filename)
    factory = javax.xml.validation.SchemaFactory.newInstance('http://www.w3.org/2001/XMLSchema');
    schema = factory.newSchema( java.io.File(schema_filename) );
    try
       schema.newValidator().validate( javax.xml.transform.stream.StreamSource( java.io.StringBufferInputStream( xmlstr ) ) );
    catch ME
       % XML string failed validation
       % Find the useful part of the error message.
       [first, last] = regexp(ME.message, '(?<=org.xml.sax.SAXParseException[:;] ).+?(?=[\r\n])');
       if ~isempty(first)
           warning('SICDSTRUCT2XML:SCHEMA_VALIDATION_FAILURE', ...
               ['XML metadata failed to validate against schema.  ' ...
               'This is not unusual.  It just means that the metadata used for ' ...
               'this dataset doesn''t exactly match the ' p1.Results.file_type ...
               ' schema being used in the writer.  In this case, ' ...
               'the (first) reason is:\n' ME.message(first:last)]);
       else % Not the error we expected
           rethrow(ME);
       end
    end
end

    %% Recursively process through metadata structure and schema to create XML
    function sicdstruct2xml_recurse(current_node, sicdmeta, schema_struct)
        if ~isempty(schema_struct)
            schema_fieldnames = fieldnames(schema_struct);
        else
            schema_fieldnames = {};
        end
        meta_fieldnames = fieldnames(sicdmeta);
        % Order fieldnames according to schema
        child_names = schema_fieldnames(ismember(schema_fieldnames, meta_fieldnames));
        % Alternative form-- only works in MATLAB 2012a and later though:
        % child_names = intersect(schema_fieldnames,meta_fieldnames,'stable');
        % Add fields in structure that are not found in schema
        child_names = [child_names; meta_fieldnames(~ismember(meta_fieldnames, child_names))];
        % Alternative 2012a+ version:
        % child_names = union(child_names,meta_fieldnames,'stable');
        % Traverse schema structure
        for i=1:numel(child_names)
            if isfield(schema_struct,child_names{i})
                child_schema_struct = schema_struct.(child_names{i});
                % Find base structure or primitive string
                while isfield(child_schema_struct, 'SCHEMA_type') && ...
                        isfield(schema_info.types, child_schema_struct.SCHEMA_type)
                    if isfield(child_schema_struct, 'SCHEMA_attributes') % Schema "extension"
                        attributes_to_pass = child_schema_struct.SCHEMA_attributes;
                    else
                        attributes_to_pass = [];
                    end
                    child_schema_struct = schema_info.types.(child_schema_struct.SCHEMA_type);
                    if ~isempty(attributes_to_pass) % Pass extension attributes to base type
                        if isfield(child_schema_struct, 'SCHEMA_attributes')
                            child_schema_struct.SCHEMA_attributes = ...
                                [child_schema_struct.SCHEMA_attributes attributes_to_pass];
                        else
                            child_schema_struct.SCHEMA_attributes = attributes_to_pass;
                        end
                    end
                end
            else
                % We try to be flexible and read all fields, regardless
                % of whether the field is described in the schema or
                % not.  This allows extra custom fields to be included that
                % may not fit the spec.  Also, if our metadata comes from
                % a SICD of a different version than the schema we are
                % using, this should allow that to work (at least
                % partially) as well.
                child_schema_struct = [];
            end

            if iscell(sicdmeta.(child_names{i})) % Cell array
                for ii = 1:numel(sicdmeta.(child_names{i}))
                    process_single_node(current_node, child_names{i}, ...
                        sicdmeta.(child_names{i}){ii}, child_schema_struct, ii);
                end
            elseif isstruct(sicdmeta.(child_names{i})) && ... % Struct array
                    numel(sicdmeta.(child_names{i}))>1
                for ii = 1:numel(sicdmeta.(child_names{i}))
                    process_single_node(current_node, child_names{i}, ...
                        sicdmeta.(child_names{i})(ii), child_schema_struct, ii);
                end
            else % Single structure/value
                child_node = process_single_node(current_node, child_names{i}, ...
                    sicdmeta.(child_names{i}), child_schema_struct);
            end

            % Add size attributes
            if isfield(child_schema_struct, 'SCHEMA_attributes')
                if any(strcmp('size', child_schema_struct.SCHEMA_attributes))
                    size_value = child_node.getChildNodes.getLength;
                    child_node.setAttribute('size',num2str(size_value));
                elseif any(strcmp('order1', child_schema_struct.SCHEMA_attributes))
                    child_node.setAttribute('order1',num2str(size(sicdmeta.(child_names{i}),1)-1));
                    if any(strcmp('order2', child_schema_struct.SCHEMA_attributes))
                        child_node.setAttribute('order2',num2str(size(sicdmeta.(child_names{i}),2)-1));
                    end
                end
            end
        end
    end % of sicdstruct2xml_recurse

    function child_node = process_single_node(current_node, node_name, sicdmeta, schema_struct, index)
        % Handle special case fieldnames first
        if any(strcmp(node_name,{'native','SICDVersion','NITF'})) % Non-spec fields added by MATLAB SAR Toolbox
            child_node = [];
        elseif strcmp(node_name,'ICP') % Special case: ICP Indexed by name rather than number
            ICP_fields = {'FRFC','FRLC','LRLC','LRFC'};
            for iii = 1:numel(ICP_fields)
                child_node = doc.createElement(node_name);
                current_node.appendChild(child_node);
                child_node.setAttribute('index',[num2str(iii) ':' ICP_fields{iii}]);
                sicdstruct2xml_recurse(child_node, sicdmeta.(ICP_fields{iii}), schema_struct)
            end
        else
            child_node = doc.createElement(node_name);
            current_node.appendChild(child_node);
            if isfield(sicdmeta, 'name') % Special attribute: Used in 'Paramater', 'Desc', 'GeoInfo', 'RefPt'
                child_node.setAttribute('name',sicdmeta.name);
                sicdmeta = rmfield(sicdmeta,'name');
            end
            if isfield(sicdmeta, 'value') % Special structure field: Used in 'Parameter', 'Desc'
                child_node.appendChild(doc.createTextNode(sicdmeta.value));
            elseif isstruct(sicdmeta) % Substructure
                sicdstruct2xml_recurse(child_node, sicdmeta, schema_struct);
                if (isfield(schema_struct,'SCHEMA_attributes') && ...
                        any(strcmp('index',schema_struct.SCHEMA_attributes))) || ...
                        (isempty(schema_struct) && exist('index','var') && ~strcmp(node_name,{'GeoInfo'}))
                    if ~exist('index','var')
                        index = 1;
                    end
                    child_node.setAttribute('index',num2str(index));
                end
            elseif (~ischar(sicdmeta) && ...% Numeric array
                    numel(sicdmeta)>1) || ...
                    (isfield(schema_struct, 'SCHEMA_attributes') && ...
                    (any(strcmp('order1', schema_struct.SCHEMA_attributes)) || ...
                    any(strcmp('size', schema_struct.SCHEMA_attributes))))
                % Either this value is passed as an array, or we know
                % it should be an array from the schema.  This test
                % catches 1x1 arrays like TimeCOAPoly for spotlight.
                order1 = size(sicdmeta,1);
                order2 = size(sicdmeta,2);
                for iii=1:order1
                    for jjj = 1:order2
                        if strcmp(node_name,'WgtFunct')
                            coef_node = doc.createElement('Wgt');
                            attribute_name = 'index';
                            val = num2str(iii);
                        elseif strcmp(node_name,'AmpTable')
                            coef_node = doc.createElement('Amplitude');
                            attribute_name = 'index';
                            val = num2str(iii-1);
                        else
                            coef_node = doc.createElement('Coef');
                            attribute_name = 'exponent1';
                            val = num2str(iii-1);
                        end
                        child_node.appendChild(coef_node);
                        coef_node.setAttribute(attribute_name,val);
                        if order2>1 || (isfield(schema_struct, 'SCHEMA_attributes') && ...
                                any(strcmp('order2', schema_struct.SCHEMA_attributes)))
                            coef_node.setAttribute('exponent2',num2str(jjj-1));
                        end
                        coef_node.appendChild(doc.createTextNode(...
                            num2str(sicdmeta(iii,jjj),'%.15E')));
                        if p1.Results.inc_class_attributes % All arrays in SICD are double type
                            child_node.setAttribute('class','xs:double');
                        end
                    end
                end
            else % Scalars
                % First check schema, then check MATLAB class of
                % value in metadata structure. If variable in
                % memory and schema disagree, we must convert type.
                if isfield(schema_struct, 'SCHEMA_type')
                    class_str = schema_struct.SCHEMA_type;
                    switch class_str
                        case 'xs:string'
                            if isnumeric(sicdmeta) % May have been incorrectly populated
                                str = num2str(sicdmeta);
                            else
                                str = sicdmeta;
                            end
                        case 'xs:double'
                            str = num2str(sicdmeta,'%.15E');
                        case {'xs:int', 'xs:integer', 'xs:positiveInteger', 'xs:nonNegativeInteger'}
                            str = num2str(sicdmeta,'%d');
                        case 'xs:dateTime'
                            if isa(sicdmeta,'double')
                                str = datestr(sicdmeta,...
                                    'yyyy-mm-ddTHH:MM:SS.FFFZ');
                            elseif isdatetime(sicdmeta)
                                % Display seconds separately since MATLAB
                                % default functions only display
                                % milliseconds.
                                str = [datestr(sicdmeta, 'yyyy-mm-ddTHH:MM') ...
                                    ':' char(regexp(num2str(sicdmeta.Second, '%012.9f'), ...
                                    ... % Trim trailing zeros except one immediately after decimal
                                    '\d*\.((\d*[^0])|(0))(?=0*$)','match')) ...
                                    'Z'];
                                % MATLAB documentation says datetime
                                % precision is good to nanoseconds so we
                                % stop at nine decimal places (11 digits if
                                % including whole seconds left of decimal).
                            end
                        case 'xs:boolean'
                            if sicdmeta
                                str = 'true';
                            else
                                str = 'false';
                            end
                        otherwise
                            error('SICDSTRUCT2XML:UNKNOWN_CLASS','Unrecognized class type.');
                    end
                else % Field not found in schema.  Guess class based on value
                    if (any(strcmp(node_name,{'DateTime','CollectStart'})) && ... % Special case: DateTime needs to be formatted/converted from double to string
                            isa(sicdmeta,'double'))
                        str = datestr(sicdmeta,...
                            'yyyy-mm-ddTHH:MM:SS.FFFZ');
                        class_str = 'xs:dateTime';
                    elseif isdatetime(sicdmeta)
                        str = [datestr(sicdmeta, 'yyyy-mm-ddTHH:MM') ...
                            ':' char(regexp(num2str(sicdmeta.Second, '%012.9f'), ...
                            ... % Trim trailing zeros except one immediately after decimal
                            '\d*\.((\d*[^0])|(0))(?=0*$)','match')) 'Z'];
                        class_str = 'xs:dateTime';
                    elseif isinteger(sicdmeta)
                        str = num2str(sicdmeta);
                        class_str = 'xs:int';
                    elseif isfloat(sicdmeta)
                        str = num2str(sicdmeta,'%.15E');
                        class_str = 'xs:double';
                    elseif islogical(sicdmeta)
                        if sicdmeta
                            str = 'true';
                        else
                            str = 'false';
                        end
                        class_str = 'xs:boolean';
                    else
                        str = sicdmeta;
                        class_str = 'xs:string';
                    end
                end
                child_node.appendChild(doc.createTextNode(str));
                if p1.Results.inc_class_attributes
                    child_node.setAttribute('class',class_str);
                end
            end
        end
    end
end % of sicdstruct2xml

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
