function s=afni_niml_print(p)
% takes a NIML data structure and converts it to string representation
%
% S=AFNI_NIML_PRINT(P) takes an NIML datastructure P and returns a string
% representation S.
%
% This function is more or less the inverse of AFNI_NIML_PARSE. 
%
% The current implementation will take any struct that has the required
% NIML fields; it does not check for any data consisentence.
%
% Please note that this function is *VERY EXPERIMENTAL*.
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>

if iscell(p)
    ss=cell(1,numel(p));
    % simple recursion
    for k=1:numel(p)
        ss{k}=afni_niml_print(p{k});
    end
    s=[ss{:}]; % concatenate
else
    headername=p.name;
    p=rmfield(p,'name'); % we don't want this as a field
    
    if isfield(p,'nodes') 
        % is a group, do recursion
        sbody=afni_niml_print(p.nodes);
        p=rmfield(p,'nodes'); % 'nodes' is not a standard NIML field
    elseif isfield(p,'data') 
        if ~isfield(p,'vec_typ')
            % in case the type is not specified, we try 
            % to derive it based on the type of p.data
            ps=derive_vec_type(p.data);
            
            % set field names if not set
            fns=fieldnames(ps);
            for k=1:numel(fns)
                if ~isfield(p,fns{k})
                    p.(fns{k})=ps.(fns{k});
                end
            end
        elseif p.vec_typ<0;
            error('vec_typ=%d not supported (yet)', p.vec_typ);
        end
                
        sbody=afni_niml_print_body(p);
        
        % some fields are not standard NIML (I think), we remove these 
        removefields=strvcat('vec_typ','vec_len','vec_num','name','data');
        fns=fieldnames(p);
        for k=1:numel(fns)
            fn=fns{k};
            if ~isempty(strmatch(fn,removefields,'exact'))
                p=rmfield(p,fn);
            end
        end
    else
        disp(p)
        error('Do not understand this struct');
    end
    
    headertext=afni_niml_print_header(p);
    s=sprintf('<%s\n%s >\n%s</%s>\n',headername,headertext,sbody,headername);
end
    
function s=afni_niml_print_body(p)

format=get_print_format(p.vec_typ,p.data);
if strcmp(format,'%s')
    around='"'; %surround by quotes if it is a string - CHECKME is that according to the standard?
else
    around='';
    p.data=p.data'; %transpose to fix major row vs. major column order
end
s=sprintf([format '\n'],p.data);

s=[around s(1:(end-1)) around]; % remove last newline
        

function s=afni_niml_print_header(p)
pre='  ';           % a bit of indendationn
post=sprintf('\n'); % every header gets a new line
fns=fieldnames(p);
n=numel(fns);

ss=cell(1,n);
for k=1:n
    fn=fns{k};
    val=strtrim(p.(fn));
    
    if isnumeric(val)
        warning('Converting numeric to string for field %s\n', fn);
        val=num2str(val);
    end
    
    % surround by quotes, if that's not the case yet
    if val(1) ~= '"' && val(end) ~= '"'
        val=['"' val '"'];
    end
    if k==n
        post=''; %no new line at the very end
    end
    ss{k}=sprintf('%s%s=%s%s',pre,fn,val,post);
end

s=[ss{:}]; % concatenate results

function p=derive_vec_type(data)
% sets the vector type, in case this is not given.
%
% data should be either a string or numeric
% TODO: support structs and cells, and mixed data types
% TODO: allow other fields missing (maybe use the whole struct rather than
%       just the data so that we can be 'smart' and use converging
%       evidence?
    nidefs=afni_ni_defs();
    if islogical(data)
        data=single(data); % convert to numeric
    end
    if isnumeric(data)
        if isequal(round(data),data)
            tp=nidefs.NI_INT;
        else
            tp=nidefs.NI_FLOAT;
        end
        [ln,nm]=size(data);
    elseif ischar(data)
        tp=nidefs.NI_STRING;
        ln=1;
        nm=1;
    else
        disp(data)
        error('Unknown data type');
    end
    
    if nm>1
        prefix=sprintf('%d*',nm);
    else
        prefix='';
    end
    
    tps=strtrim(nidefs.type_name(tp+1,:)); % tricky base0 vs base1
    
    p.ni_type=[prefix tps];
    p.ni_dimen=sprintf('%d',ln); %NNO Jan 2010 fix
    
    p.vec_typ=repmat(tp,1,nm); % tricky base0 vs base1 % Jan 2010: removed 'tp+1'
    p.vec_len=ln;
    p.vec_num=nm;


function f=get_print_format(vec_typ,data,nidefs)
% use vec_typ to create a format string for printf
% data is used in case we print floats; we use the least possible number of
% floating point positions

    if nargin<3
        nidefs=afni_ni_defs();
    end

    n=numel(vec_typ);
    if n>1
        if size(data,2) ~= n
            error('Mismatch between data and vec_typ: %d ~= %d', n, size(data,2));
        end
        fs=cell(1,n);
        for k=1:n
            if k<n
                addend=' '; % spacing between elements ...
            else
                addend='';  % ... except for the last element
            end
            % recursion
            fs{k}=[get_print_format(vec_typ(k),data(:,k),nidefs) addend];
        end
        f=[fs{:}]; % concatenate result
    else
        vec_str=strtrim(nidefs.type_name(vec_typ+1,:)); %tricky again, base0 vs base1

        switch vec_str
            case 'String'
                f='%s';
            case {'byte','short','int'} % CHECKME don't know if this is ok for bytes and shorts...
                f='%d';
            otherwise
                precision=get_required_precision(data);
                f=sprintf('%%.%df',precision);
        end
    end
    

        
    