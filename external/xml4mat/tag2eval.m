function [to_eval,eval_i,j]=tag2eval(tag_contents,to_eval,eval_i,j,tags)

%TAG2EVAL Extacts statements for evaluation by XML2MAT
%
%Syntax: [to_eval,eval_i,j]=tag2eval(tag_contents,to_eval,j)
%
%Description:
% Autorecursive function that parses tag_contents structure
% generated within XML2MAT (see %RECOVER STRUCTURE while loop)
%
% See Also: xml2whos
%
% Jonas Almeida, almeidaj@musc.edu 20 Aug 2002, XML4MAT Tbox

i=tags(j);
w=tag_contents{i,2};%disp(i)
if ~isfield(w,'class');w.class='struct';w.size=[1 1];tag_contents{i,2}=w;end
if strcmp(w.class,'struct')
    if ~isfield(w,'size')
        n_fields=length(w.fields);
        for f_i=1:n_fields
            field_names{f_i}=tag_contents{w.fields(f_i),2}.name;
        end
        unique_names=unique(field_names);
        n_unique=length(unique_names);
        n_certo=(n_fields/n_unique);
        w.size=[1 n_certo];
    end
    nn=prod(w.size); %number of elements
    nf=length(w.fields); %number of fields per element
    iis='i1';for ii=2:length(w.size);iis=[iis,',i',num2str(ii)];end %indexes
    eval(['[',iis,']=ind2sub(w.size,[1:nn]);']); % assigning values to indexes
    iis='i1(ind)';for ii=2:length(w.size);iis=[iis,',i',num2str(ii),'(ind)'];end % indexes of indexes
    %disp(iis)
    for ind=1:nn    
        for ind_f=1:(nf/nn)
            %disp(['ind=',num2str(ind),' nf=',num2str(nf)])
            j=j+1;%disp(['j=',num2str(j)])
            i=tags(j);
            field_name=tag_contents{w.fields(ind_f),2}.name;iis_val=num2str(eval(['[',iis,']']));
            iis_val(findstr(iis_val,'  '))=[];iis_val(isspace(iis_val))=',';
            eval_i{length(eval_i)+1}=['(',iis_val,').',field_name];
            [to_eval,eval_i,j]=tag2eval(tag_contents,to_eval,eval_i,j,tags);
            %eval_i_str='';for eval_i_j=1:length(eval_i);eval_i_str=[eval_i_str,eval_i{eval_i_j}];end;eval_i_str=[eval_i_str,'=tag_contents{',num2str(i),',4};'];disp(eval_i_str)            
            %to_eval{length(to_eval)+1}=eval_i_str;
            eval_i(end)=[];
        end
    end
elseif strcmp(w.class,'cell')
    if ~isfield(w,'size')
        w.size=[1,length(w.fields)];
    end
    nn=prod(w.size); %number of elements
    nf=length(w.fields); %number of fields per element
    iis='i1';for ii=2:length(w.size);iis=[iis,',i',num2str(ii)];end %indexes
    eval(['[',iis,']=ind2sub(w.size,[1:nn]);']); % assigning values to indexes
    iis='i1(ind)';for ii=2:length(w.size);iis=[iis,',i',num2str(ii),'(ind)'];end % indexes of indexes
    %disp(iis)
    for ind=1:nn
        for ind_f=1:(nf/nn)
            %disp(['ind=',num2str(ind),' nf=',num2str(nf)])
            j=j+1;%disp(['j=',num2str(j)])
            i=tags(j);
            %disp(iis)
            
            field_name=tag_contents{w.fields(ind_f),2}.name;iis_val=num2str(eval(['[',iis,']']));
            iis_val(findstr(iis_val,'  '))=[];iis_val(isspace(iis_val))=',';
            eval_i{length(eval_i)+1}=['{',iis_val,'}'];
            [to_eval,eval_i,j]=tag2eval(tag_contents,to_eval,eval_i,j,tags);
            %eval_i_str='';for eval_i_j=1:length(eval_i);eval_i_str=[eval_i_str,eval_i{eval_i_j}];end;eval_i_str=[eval_i_str,'=tag_contents{',num2str(i),',4};'];disp(eval_i_str)            
            %to_eval{length(to_eval)+1}=eval_i_str;
            eval_i(end)=[];
        end
    end
elseif strcmp(w.class,'cellstruct')
    %disp(w)
    eval_i{end+1}={};
    for i=1:length(w.fields)
        ww=tag_contents{w.fields(i),2};
        field_exist{i}=ww.name;
        n=sum(strcmp(ww.name,field_exist));
        %disp(['field [',ww.name,'] ',num2str(n)]);
        %eval_i{length(eval_i)+(i==1)}=['.',ww.name,'{',num2str(n),'}'];
        eval_i{length(eval_i)}=['.',ww.name,'{',num2str(n),'}'];
        j=j+1;
        [to_eval,eval_i,j]=tag2eval(tag_contents,to_eval,eval_i,j,tags);
        %if ((i>1)&(n==1));to_eval(end)=[];end
        %disp([eval_i,'.',ww.name])
        %disp(BM_tag2eval(
    end
    eval_i(end)=[];
    clear field_exist;
else % str or double
    %disp('not struct')
    eval_i_str='';for eval_i_j=1:length(eval_i);eval_i_str=[eval_i_str,eval_i{eval_i_j}];end;eval_i_str=[eval_i_str,'=tag_contents{',num2str(i),',4};'];%disp(eval_i_str)
    to_eval{length(to_eval)+1}=eval_i_str;
end
    



%j=j+1;