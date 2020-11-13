function CCF = parseCCF(varargin)

% parseCCF
%
% Parses an XML CCF file.
%
%   Nick Halper
%   support@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.1.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
% 1.1.0.0: January 18, 2016 - Kian Torab
%   - Minor bug fix with file loading.
%
% 1.1.1.0: January 19, 2016 - Kian Torab
%   - Added a progress bar.
%
% 1.1.2.0: October 20, 2016 - Saman Hagh-gooie
%   - Fixed a invalid character bug.
%   - Bug fixes with file loading
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 0
    fullfilename = varargin{1};
else
    fullfilename = [];
end

if ~exist(fullfilename, 'file')
    [fileName pathName] = getFile('*.ccf', 'Choose a CCF file...');
    fullfilename = [pathName fileName];
else
    [pathName, fileName, ext] = fileparts(fullfilename);
    pathName = [pathName];
    fileName = [fileName, ext];
end

prevDecodedXMLCCF = [fullfilename(1:end-3) 'xmld'];

if exist(prevDecodedXMLCCF, 'file') == 2
    load(prevDecodedXMLCCF, '-mat');
    return;
end

disp('The CCF loading can take up to half an hour the first time. Please be patient.');

OBJ = xmlread(fullfile(pathName,fileName));
removeIndentNodes(OBJ.getChildNodes);
CCF = parseChildNodes(OBJ);

save(prevDecodedXMLCCF, 'CCF');

%function children = parseChildNodes(OBJ)
function children = parseChildNodes(theNode)
% Recurse over node children.
children = 0;
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);    

   counter=0; % added by SH 05.oct.2016
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
        counter = counter + 1;
        if mod(counter,20) == 0
            fprintf('.');
        end
    end
end


%function nodeStruct - makeStructFromNode(OBJ)
function nodeStruct = makeStructFromNode(theNode)
nodeStruct = struct(...
    'Name',char(theNode.getNodeName),...
    'Attributes', parseAttributes(theNode),...
    'Data','',...
    'Children', parseChildNodes(theNode));


if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end




function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end


function removeIndentNodes( childNodes )

numNodes = childNodes.getLength;
remList = [];
counter = 0;
for i = numNodes:-1:1
    counter = counter + 1;
    if rem(counter,20) == 0
        fprintf('.');
    end
	theChild = childNodes.item(i-1);
	if (theChild.hasChildNodes)
        removeIndentNodes(theChild.getChildNodes);
    else
        if ( theChild.getNodeType == theChild.TEXT_NODE && ...
           ~isempty(char(theChild.getData()))         && ...
           all(isspace(char(theChild.getData()))))
         remList(end+1) = i-1; % java indexing
        end
    end
end
for i = 1:length(remList)
   childNodes.removeChild(childNodes.item(remList(i)));
end



