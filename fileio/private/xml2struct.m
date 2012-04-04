function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</DifferentElement>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Used to produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Will produce (gp: to matche the output of xml2struct in XML4MAT, but note that Element(2) is empty):
% Element: Some text 
% DifferentElement:
%    attrib2: 2 
%    DifferentElement: Some more text 
% attrib1: Some value 
% 
% Element:  
% DifferentElement:
%    attrib3: 2 
%    attrib4: 1 
%    DifferentElement: Even more text 
% attrib1: 
%
% Note the characters : - and . are not supported in structure fieldnames and
% are replaced by _
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% 2011/12/14 giopia: changes in the main function to make more similar to xml2struct of the XML4MAT toolbox, bc it's used by fieldtrip
% 2012/04/04 roboos: added the original license clause, see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=645#c11
% 2012/04/04 roboos: don't print the filename that is being read

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2010, Wouter Falkena
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (nargin < 1)
        clc;
        help xml2struct
        return
    end
    
    %check for existance
    if (exist(file,'file') == 0)
        %Perhaps the xml extension was omitted from the file name. Add the
        %extension and try again.
        if (isempty(strfind(file,'.xml')))
            file = [file '.xml'];
        end
        
        if (exist(file,'file') == 0)
            error(['The file ' file ' could not be found']);
        end
    end
    %fprintf('xml2struct reading %s\n', file); % gp 11/12/15
    %read the xml file
    xDoc = xmlread(file);
    
    %parse xDoc into a MATLAB structure
    s = parseChildNodes(xDoc);
    fn = fieldnames(s); s = s.(fn{1}); % gp 11/12/15: output is compatible with xml2struct of xml2mat
end

% ----- Subfunction parseChildNodes -----
function [children,ptext] = parseChildNodes(theNode)
    % Recurse over node children.
    children = struct;
    ptext = [];
    if theNode.hasChildNodes
        childNodes = theNode.getChildNodes;
        numChildNodes = childNodes.getLength;

        for count = 1:numChildNodes
            theChild = childNodes.item(count-1);
            [text,name,attr,childs] = getNodeData(theChild);

            if (~strcmp(name,'#text') && ~strcmp(name,'#comment'))
                %XML allows the same elements to be defined multiple times,
                %put each in a different cell
                if (isfield(children,name))
%                     if 0 % numel(children) > 1 % gp 11/12/15: (~iscell(children.(name)))
%                         %put existsing element into cell format
%                         children.(name) = {children.(name)};
%                     end
                    index = length(children)+1; % gp 11/12/15: index = length(children.(name))+1;
                else 
                  index = 1; % gp 11/12/15: new field
                end
                %add new element
                children(index).(name) = childs; 
                if isempty(attr)
                  if(~isempty(text))
                    children(index).(name) = text; 
                  end
                else
                  fn = fieldnames(attr);
                  for f = 1:numel(fn)
                    children(index).(name)(1).(fn{f}) = attr.(fn{f}); % gp 11/12/15: children.(name){index}.('Attributes') = attr;
                  end
                  if(~isempty(text))
                    children(index).(name).(name) = text; % gp 11/12/15: children.(name){index}.('Text') = text;
                  end
                end

%                 else % gp 11/12/15: cleaner code, don't reuse the same code
%                     %add previously unknown new element to the structure
%                     children.(name) = childs;
%                     if(~isempty(text)) 
%                         children.(name) = text; % gp 11/12/15: children.(name).('Text') = text; 
%                     end
%                     if(~isempty(attr)) 
%                         children.('Attributes') = attr; % gp 11/12/15 children.(name).('Attributes') = attr; 
%                     end
%                 end
            elseif (strcmp(name,'#text'))
                %this is the text in an element (i.e. the parentNode) 
                if (~isempty(regexprep(text,'[\s]*','')))
                    if (isempty(ptext))
                        ptext = text;
                    else
                        %what to do when element data is as follows:
                        %<element>Text <!--Comment--> More text</element>
                        
                        %put the text in different cells:
                        % if (~iscell(ptext)) ptext = {ptext}; end
                        % ptext{length(ptext)+1} = text;
                        
                        %just append the text
                        ptext = [ptext text];
                    end
                end
            end
        end
    end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs] = getNodeData(theNode)
    % Create structure of node info.
    
    %make sure name is allowed as structure name
    name = regexprep(char(theNode.getNodeName),'[-:.]','_');

    attr = parseAttributes(theNode);
    if (isempty(fieldnames(attr))) 
        attr = []; 
    end
    
    %parse child nodes
    [childs,text] = parseChildNodes(theNode);
    
    if (isempty(fieldnames(childs)))
        %get the data of any childless nodes
        try
            %faster then if any(strcmp(methods(theNode), 'getData'))
            text = char(theNode.getData);
        catch
            %no data
        end
    end
    
end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.

    attributes = struct;
    if theNode.hasAttributes
       theAttributes = theNode.getAttributes;
       numAttributes = theAttributes.getLength;

       for count = 1:numAttributes
            %attrib = theAttributes.item(count-1);
            %attr_name = regexprep(char(attrib.getName),'[-:.]','_');
            %attributes.(attr_name) = char(attrib.getValue);

            %Suggestion of Adrian Wanner
            str = theAttributes.item(count-1).toString.toCharArray()'; 
            k = strfind(str,'='); 
            attr_name = regexprep(str(1:(k(1)-1)),'[-:.]','_'); 
            attributes.(attr_name) = str((k(1)+2):(end-1));
       end
    end
end
