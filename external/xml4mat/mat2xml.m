function XML=mat2xml(MAT,VARNAME)

% MAT2XML converts structured variable MAT into XML string
%
%Syntax XML=mat2xml(MAT,VARNAME)
%
%Description 
%  MAT : structured varable
%  VARNAME : variable name (string)
%  XML : xml version of structured variable (string)
%
% See Also: XML2MAT
%
% Jonas Almeida, almeidaj@musc.edu, 20 Aug 2002, XML4MAT Tbox

if nargin<2;VARNAME='ans';end % if not provided make it a matlab answer variable
w=whos('MAT');w.name=VARNAME;

XML=['<',w.name,' class="',w.class,'" size="',num2str(w.size),'">'];
if strcmp(w.class,'char')
    XML=[XML,spcharin(MAT(:)')];
elseif strcmp(w.class,'struct')
    names=fieldnames(MAT);
    %struct_fields=[' fields="',names{1}];for j=2:length(names);struct_fields=[struct_fields,' ',names{j}];end;struct_fields=[struct_fields,'">'];XML=[XML(1:(end-1)),struct_fields];
    for i=1:prod(w.size)
        for j=1:length(names)
            XML=[XML,mat2xml(eval(['MAT(i).',names{j}]),names{j})];
        end
    end
elseif strcmp(w.class,'cell')
    for i=1:prod(w.size)
        XML=[XML,mat2xml(MAT{i},'cell')];        
    end
else %if strcmp(w.class,'double')|strcmp(w.class,'single')
    XML=[XML,num2str(MAT(:)')];    
end
XML=[XML,'</',w.name,'>'];
XML=regexprep(XML,'\s{2,}',' ');
