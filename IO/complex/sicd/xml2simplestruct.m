function [ output_XML_struct ] = xml2simplestruct( input_domnode )
%XML2SIMPLESTRUCT Convert an XML Document Object Model (DOM) to a simple MATLAB structure
%
% This structure contains only the text element tree of the XML with its
% attibutes.  This won't contain all of the power and flexibility of XML,
% but it is easier to browse and reference in MATLAB for many simple XML
% structures.
%
% This routine doesn't handle CDATA, processing instructions, comments, or
% anything more complicated that a tree of named text elements with
% attributes.  We also hope that there are no tags named "attributes" in
% our XML file.  But for most basic XML files, this should work fine.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% XML must have single root (document node).  Start below that for real
% content
output_XML_struct=simplify_node(input_domnode.getFirstChild);

    function out_node = simplify_node( in_node )
        out_node = struct;
        if in_node.hasChildNodes
            childNodes = in_node.getChildNodes;
            numChildNodes = childNodes.getLength;
            for i=1:numChildNodes
                current_child=childNodes.item(i-1);
                current_name=genvarname(char(current_child.getNodeName));
                if current_child.getNodeType==3 % Text node.  Save value
                   % Save only text values.  Discard other data types (CDATA, comments, etc.)
                   out_node.data = char(current_child.getData);
                elseif current_child.getNodeType==1 % Element node
                    if ~isfield(out_node,current_name) % First element with this tag
                        out_node.(current_name)=simplify_node(current_child);
                    else % Multiple elements with same tag
                        if ~iscell(out_node.(current_name)) % Convert to cell array if more than one
                            out_node.(current_name)={out_node.(current_name)};
                        end
                        out_node.(current_name){end+1}=simplify_node(current_child); % Append to end
                    end
                end
            end
        end
        if in_node.hasAttributes % Add list of attributes in their own field
           attributeList = in_node.getAttributes;
           for i = 1:attributeList.getLength
              current_attribute = attributeList.item(i-1);
              out_node.attributes.(genvarname(char(current_attribute.getName)))=char(current_attribute.getValue);
           end
        end
    end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////