function xmlStruct = module_read_neurone_xml(xmlfile)
%MODULE_READ_NEURONE_XML   Read NeurOne XML-files into structures.
%
%  Input  : A NeurOne XML file (e.g. Session.xml, Protocol.xml,
%           Triggers.xml)
%  Output : A structure with the information from the XML-file.
%  Example: Session  = read_neurone_xml('/data/Session.xml');
%           Protocol = read_neurone_xml('/data/Protocol.xml');
%           Triggers = read_neurone_xml('/data/Triggers.xml');
%
%  Dependencies: none
%
%  Module_read_neurone_xml is part of NeurOne Tools for Matlab.
%
%  The NeurOne Tools for Matlab consists of the functions:
%       module_read_neurone.m, module_read_neurone_data.m,
%       module_read_neurone_events.m, module_read_neurone_xml.m
%
%  ========================================================================
%  COPYRIGHT NOTICE
%  ========================================================================
%  Copyright 2009, 2010 Andreas Henelius (andreas.henelius@ttl.fi)
%  Finnish Institute of Occupational Health (http://www.ttl.fi/)
%  ========================================================================
%  This file is part of NeurOne Tools for Matlab.
% 
%  NeurOne Tools for Matlab is free software: you can redistribute it
%  and/or modify it under the terms of the GNU General Public License as
%  published by the Free Software Foundation, either version 3 of the
%  License, or (at your option) any later version.
%
%  NeurOne Tools for Matlab is distributed in the hope that it will be
%  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with NeurOne Tools for Matlab.
%  If not, see <http://www.gnu.org/licenses/>.
%  =======================================================================

%% Read XML-data using the builtin xmlread-function and get all child nodes
% of the root XML-element.
xml = xmlread(xmlfile);
children = xml.getChildNodes;

%% Begin walking through the XML-tree
% We initialise the walk with the children of the xml root element, set level
% at zero and give an empty structure to be filled with data from the xml.
xmlStruct = walkNodes(children, 0, {});

% Recursive function for walking through xml document.
    function xmlStruct = walkNodes(node, level, xmlStruct)
        if level < 3
            if (node.hasChildNodes)
                moreChildren = node.getChildNodes;
                for i=0:moreChildren.getLength-1
                    % 'still has child nodes, going on'
                    nextNode = moreChildren.item(i);
                    nodeName = char(nextNode.getNodeName);
                    % Skip "empty" nodes (#text)
                    if ~strcmpi(nodeName,'#text')
                        if nextNode.hasChildNodes
                            % Read node data
                            nodeData = node2struct(nextNode);
                            if ~isempty(fieldnames(nodeData))
                                % If we are at recursion depth 1, use the
                                % nodeNames for main field names of the
                                % structure to be returned. Append data to the
                                % xml-structure and create subfields if the
                                % field already exists.
                                if level==1
                                    if ~isfield(xmlStruct, nodeName)
                                        xmlStruct(1).(nodeName) = nodeData;
                                    else
                                        xmlStruct.(nodeName)(end+1) = nodeData;
                                    end
                                end
                            end
                        end
                    end
                    % Some childnodes left, recurse to access them.
                    xmlStruct = walkNodes(moreChildren.item(i), level + 1, xmlStruct);
                end
            else
                % 'no more child nodes'
            end
        end
    end


%% Helper function to return data from a node
    function nodeData = node2struct(mainNode)
        nodeData = struct();
        childNodes = mainNode.getChildNodes;
        for n=0:childNodes.getLength-1
            childNode = childNodes.item(n) ;
            nodeName = char(childNode.getNodeName);
            % Skip "empty" nodes (#text)
            if ~strcmpi(nodeName,'#text') && isempty(strfind(nodeName,':'))
                % If the node has children (i.e. the element has a value),
                % return possible data, otherwise return empty string
                % if the element has no value.
                if childNode.hasChildNodes
                    nodeData.(nodeName) = char(childNode.getFirstChild.getNodeValue);
                else
                    nodeData.(nodeName) = '';
                end
            end
        end
    end
end % end of read_neurone_xml.m
