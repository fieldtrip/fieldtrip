function [MAT,VARNAME,tag_contents]=xml2mat(XML)

%XML2MAT converts XML string into matlab structure variable
%Syntax: [MAT,VARNAME,tag_contents]=xml2mat(XML)
%Description
% XML is an XML formated string using rules compliant with
%     the proceedure implemented in MAT2XML.
%     It can also be a file name with MbML text.
% VARNAME is the variable name.
%     If recovering the MAT variable with the original name is desired than
%     the followinf line will do the trick:
%
%          [MAT,VARNAME]=xml2mat(XML);eval([VARNAME,'=MAT'])
%
% See Also: mat2xml
% Jonas Almeida, 20 Aug 2002, XML4MAT Tbox


if strncmp(XML,'%60;',4)  % XML provided as encoded XML string
    XML=spcharout(XML);
end

% is this a xml string or xml filename ?
if XML(1)~='<' 
    XML=strrep(file2str(XML),'''','''''');
    % Remove non-content lines if they exist
    XML=regexprep(XML,'<[?!].*?>','');
end

%Analise XML line
tag_ini=find(XML=='<');
tag_end=find(XML=='>');
n=length(tag_ini); % number of enclosed tag structures

% extract tag_names properties and contents
if n>0
for i=1:n
    tag_close(i)=(XML(tag_ini(i)+1)=='/'); % 1 for closing contents and 0 for opening
    tag_contents{i,1}=XML(tag_ini(i)+1+tag_close(i):tag_end(i)-1);  % first column contains names
end
tag_path=[0]; % first name is root name
to_do=zeros(n,1); % 1 needs doing
to_do(1)=1;
tag_open=~tag_close;
i=1;tag_contents{i,2}=xml2whos(tag_contents{i,1});tag_contents{i,2}.fields=[];tag_contents{i,3}=tag_path;tag_path=[tag_path,i];
for i=2:n
    tag_contents{i,2}=xml2whos(tag_contents{i,1});tag_contents{i,2}.fields=[]; % second column contains WHO properties
    if tag_open(i)
        tag_contents{i,3}=tag_path;  % third column contains structre path
        tag_path=[tag_path,i];
        to_do(i)=1; % do this one
        to_do(tag_contents{i,3}+(tag_contents{i,3}==0))=0; % do host later
        tag_contents{tag_path(end-1),2}.fields=[tag_contents{tag_path(end-1),2}.fields,i];
    else
        tag_path(end)=[]; % move back one level
    end  
end

% RECOVER DATA
todo_list=find(to_do==1)';do_i=1;
while do_i<=length(todo_list);
    i=todo_list(do_i);
    % for each case extract value
    tag_contents{i,4}=XML(tag_end(i)+1:tag_ini(i+1)-1);  % 4th column contains value as string
    w=tag_contents{i,2};
    % recover number format if appropriated
    if ~isfield(w,'class')
        w.class='char';
        tag_contents{i,4}=spcharout(tag_contents{i,4});
        w.size=size(tag_contents{i,4}); %correct size tag as well
        %tag_contents{i,2}=w;
    elseif strcmp(w.class,'char')
        tag_contents{i,4}=spcharout(tag_contents{i,4});
        w.size=size(tag_contents{i,4}); %correct size tag as well
        %tag_contents{i,2}=w;
    elseif strcmp(w.class,'struct')&(~isfield(w,'size'))
        n_fields=length(w.fields);
        for i=1:n_fields
            field_names{i}=tag_contents{w.fields(i),2}.name;
        end
        unique_names=unique(field_names);
        n_unique=length(unique_names);
        n_certo=(n_fields/n_unique);
        w.size=[1 n_certo];
    else % it is a numeric type, say "double" or "single"
        tag_contents{i,4}=str2num(tag_contents{i,4});
        if ~isfield(w,'size');
            w.size=size(tag_contents{i,4});
        end
    end
    tag_contents{i,2}=w;
    if (length(w.size>2)|(w.size(1)>1))
        iis='i1';for ii=2:length(w.size);iis=[iis,',i',num2str(ii)];end
        nn=prod(w.size); %number of elements
        eval(['[',iis,']=ind2sub(w.size,[1:nn]);']); % generation of indexes
        iis='i1(ind)';for ii=2:length(w.size);iis=[iis,',i',num2str(ii),'(ind)'];end % indexes of indexes
        for ind=1:nn
            eval(['valor(',iis,')=tag_contents{i,4}(ind);'])            
        end        
        if exist('valor')==1;tag_contents{i,4}=valor;clear valor;end
    end
    do_i=do_i+1;
end

% RECOVER STRUCTURE
j=1;tags=find(tag_open);
eval_i={'MAT'};to_eval={};
while j<=length(tags)
    i=tags(j);
    [to_eval,eval_i,j]=tag2eval(tag_contents,to_eval,eval_i,j,tags);
    %w=tag_contents{i,2}    
    j=j+1;
end

for i=1:length(to_eval)
    %disp(to_eval{i})
    eval(to_eval{i})
end

%MAT.tag_contents=tag_contents;
VARNAME=tag_contents{1,2}.name;
else
    MAT=XML;
    VARNAME=[];
    tag_contents=[];
end