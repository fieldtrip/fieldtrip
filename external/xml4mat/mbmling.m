function yml=mbmling(xml,plt)

%MBMLING_CELL converts any XML into MbML (Matlab Markup Language) compliant string
%
%Syntax: yml=mbmling(xml,plt)
%
%Description:
%   Converts any XML syntax into Matlab Markup Language (MbML)
%
% Jonas Almeida, almeidaj@musc.edu, 18 May 2002, MAT4NAT Tbox

if nargin<2;plt=0;end
if plt==1;disp('MbMLing progress report');disp('---------BEGIN----(9 steps)----------');end

% Remove non-content lines if they exist
if plt==1;disp('1. Removing non-content lines if they exist');end
xml=regexprep(xml,'<[?!].*?>','');
% Remove empty spaces between tags
if plt==1;disp('2. Removing empty spaces between tags');end
%xml=regexprep(xml,'>[ ]+?<','><');
xml=strrep(xml,char([13,10]),' '); %replace chariage returns by spaces
n=length(xml)+1;
while n>length(xml)
    n=length(xml);
    xml=strrep(xml,'> ','>');
end
% Replace symbols that may conflict with matlab variable naming by underscore characters
if plt==1;disp('3. Replacing symbols that may conflict with matlab variable naming by underscore characters');end
xml=regexprep(xml,'<[\/]{0,1}(\w*)[^\w\/> ]+(.*?)([ >])','<$1_$2$3','tokenize');
% Replace one-tag format by open / close tagging
if plt==1;disp('4. Replace one-tag format by open / close tagging');end
xml_=regexprep(xml,'<(\w+)([^>]*?)/>','<$1$2></$1>','tokenize');
% Turn attributes into contents
if plt==1;disp('5. Turn attributes into contents');end
xml=regexprep(xml_,'<([^>]+) +(\w+)="(.*?)" *>','<$1><$2>$3</$2>','tokenize');
while ~strcmp(xml,xml_)
    %disp('...')
    xml_=xml;xml=regexprep(xml_,'<([^>]+) +(\w+)="(.*?)" *>','<$1><$2>$3</$2>','tokenize');
end
% remove leftover spaces in tag names
if plt==1;disp('6. Removing leftover spaces in tag names');end
xml=regexprep(xml,'(<\w+) +>','$1>','tokenize');
% Tag untagged contents
if plt==1;disp('7. Tagging untagged contents');end
yml=regexprep(xml,'(</\w+>)([^<>]+)</(\w+)','$1<$3>$2</$3></$3','tokenize');
% Tag cell structures
if plt==1;disp('8. Tag cell structures');end
yml=regexprep(yml,'<(\w+)>','<$1 class="cell"><cell>','tokenize');
yml=regexprep(yml,'<(/\w+)>','</""""><$1>','tokenize');
yml=regexprep(yml,'<(/\w+)>(<\w+ )','<$1></""""><cell>$2','tokenize');
% Remove cell tag arround content values
if plt==1;disp('9. Remove cell tag arround content values');end
yml=regexprep(yml,'class="cell"><cell>([^<]*)</"""">','class="char">$1','tokenize');
yml=strrep(yml,'</"""">','</cell>');
if plt==1;disp('------------END----------------------');end
