function w=xml2whos(w_xml)

%XML2WHOS identifies WHOS-type structured variable from xml descriptor
%
%Syntax: w=xml2whos(w_xml)
%
%Description
% w is a structured variable with the information that would have benn
%   returned by a WHOS command of the structured variable
% w_xml is the corresponding xml descriptor, e.g. '<name class size>'
%
% See Also: xml2mat, whos
%
% Jonas Almeida, almeidaj@mussc.edu, 20 Aug 2002, XML4MAT Tbox

w_xml=[w_xml,' '];
% extract name
i=1;while ~isspace(w_xml(i));i=i+1;end;w.name=w_xml(1:i-1);
% extract properties
i=i+1;i_1=i;is_value=0;
while i<length(w_xml)
    while ((~isspace(w_xml(i)))|is_value);i=i+1;if w_xml(i)=='"';is_value=~is_value;end;end;i_end=i-1;descrp=w_xml(i_1:i_end); % extract description
    j=1;while descrp(j)~='=';j=j+1;end;propt=descrp(1:j-1); %extract property name
    j=j+1;if descrp(j)~='"';error(['XML error 03: property value is not delimited by " " :',w_xml]);end;j=j+1;
    j_1=j;while descrp(j)~='"';j=j+1;end;j_end=j-1;propt_value=descrp(j_1:j_end); %extract property value
    %disp(['descriptor: ',descrp]);disp(['name: ',propt]);disp(['value: ',propt_value])
    eval(['w.',propt,'=''',propt_value,''';'])
    i=i+1;i_1=i;
end

if isfield(w,'size');w.size=str2num(w.size);end