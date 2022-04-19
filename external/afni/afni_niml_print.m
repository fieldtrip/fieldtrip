function s=afni_niml_print(p, format)
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

    if nargin<2
        format='binary';
    end

    format=update_output_format(format);
    s=afni_niml_print_helper(p, format);


function format=update_output_format(format)
    supported_formats={'ascii','binary',...
                        'binary.lsbfirst','binary.msbfirst'};

    if strcmp(format,'binary')
        [unused,unused,endian]=computer();
        format=sprintf('binary.%ssbfirst',lower(endian));
    end

    if isempty(strmatch(format, supported_formats))
        error('illegal format %s, use one of:%s', ...
                format, sprintf(' %s', supported_formats{:}));
    end




function s=afni_niml_print_helper(p, format)
    if iscell(p)
        ss=cell(1,numel(p));
        % simple recursion
        for k=1:numel(p)
            ss{k}=afni_niml_print_helper(p{k},format);
        end
        s=[ss{:}]; % concatenate
    else
        headername=p.name;
        p=rmfield(p,'name'); % we don't want this as a field

        if isfield(p,'nodes')
            % is a group, do recursion
            sbody=uint8(afni_niml_print_helper(p.nodes,format));

            field_to_remove={'nodes'};

        elseif isfield(p,'data')

            % in case the type is not specified, we try
            % to derive it from the data
            p=set_niml_vec_typ_attributes(p);

            [sbody,ni_form]=afni_niml_print_data(p,format);
            if ~isempty(ni_form)
                p.ni_form=ni_form;
            end

            % some fields are not standard NIML (I think), we remove these
            field_to_remove={'vec_typ','vec_len','vec_num','name','data'};

        else
            error('Do not understand input from class %s', class(p));
        end

        to_remove=intersect(fieldnames(p),field_to_remove);
        p=rmfield(p, to_remove);

        headertext=afni_niml_print_header(p);

        prefix=uint8(sprintf('<%s\n%s >',headername,headertext));
        postfix=uint8(sprintf('</%s>\n',headername));

        s=[prefix sbody postfix];
    end


function [s, ni_form]=afni_niml_print_data(p,format)

    out_format_is_ascii=strcmp(format,'ascii');
    data_is_all_numeric=vec_typ_and_data_is_all_numeric(p.vec_typ, p.data);

    if out_format_is_ascii || ~data_is_all_numeric
        s=afni_niml_print_body_ascii(p);
        ni_form=[];
    else
        s=afni_niml_print_body_binary(p,format);
        ni_form=format;
    end


function string_data=afni_niml_print_body_ascii(p)
    printer=get_ascii_printer(p.vec_typ,p.data);
    string_data=uint8(printer(p.data));




function binary_data=afni_niml_print_body_binary(p, format)
    vec_typ=p.vec_typ;
    if numel(vec_typ)>1 && any(vec_typ(1)~=vec_typ)
        error('binary only supported for uniform vec_typ');
    end

    ni_defs=afni_ni_defs();
    converter=ni_defs.type{vec_typ(1)+1}; % base 1
    data=converter(p.data');

    [unused,unused,endian]=computer();
    computer_format=sprintf('binary.%ssbfirst',lower(endian));

    do_swap=~strcmp(format, computer_format);

    if do_swap
        data=afni_swapbytes(data);
    end

    binary_data=typecast(data(:)','uint8');




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
            val=num2str(val);
        elseif ~ischar(val)
            error('expected char or numeric value, found %s', class(val));
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

function niml=set_niml_vec_typ_attributes(niml)
% sets the vector type, in case this is not given.
    if ~isfield(niml,'vec_typ') || ~isfield(niml,'vec_len')
        [niml.vec_typ, niml.vec_len]=get_vec_typ_len_from_data(...
                                                            niml.data);
    end
    niml.vec_num=numel(niml.vec_typ);

    niml.ni_type=vec_typ2string(niml.vec_typ);
    niml.ni_dimen=sprintf('%d',niml.vec_len);

function ni_type=vec_typ2string(vec_typ)
    % converts numeric vec_typ to string representation
    ni_defs=afni_ni_defs();
    type_name=cellstr(ni_defs.type_name);

    n=numel(vec_typ);
    ni_type_cell=cell(1,2*n);

    pos=0;

    count=0;
    for k=1:n
        tp=vec_typ(k);
        count=count+1;
        next_is_same=k<n && vec_typ(k+1)==tp;

        if next_is_same
            continue;
        end

        assert(count>0);

        if count>1
            pos=pos+1;
            ni_type_cell{pos}=sprintf('%d*',count);
        end

        pos=pos+1;
        ni_type_cell{pos}=type_name{tp+1};

        count=0;
    end

    ni_type=cat(2,ni_type_cell{1:pos});


function [vec_typ, vec_len]=get_vec_typ_len_from_data(data)
    % helper function to get vec type from data
    if iscell(data) && ~iscellstr(data)
        % do recursive call
        vec_typ=cellfun(@get_vec_typ_len_from_data,data);

        if ~isvector(vec_typ)
            error('Unrecognized data');
        end

        vec_len_s=cellfun(@numel,data);
        if ~all(vec_len_s(1)==vec_len_s)
            error('data has different sizes');
        end

        vec_len=vec_len_s(1);

        return
    end

    [vec_len,vec_num]=size(data);

    if islogical(data)
        data=single(data); % convert to numeric
    end

    ni_defs=afni_ni_defs();

    type_string_s=cellfun(@func2str, ni_defs.type, 'UniformOutput',false);



    tp=[];


    if isnumeric(data)
        if isequal(round(data),data)
            tp=ni_defs.NI_INT;
        else
            tp=ni_defs.NI_FLOAT;
        end

    elseif ischar(data) || iscellstr(data)
        tp=ni_defs.NI_STRING;

    else
         % check builtin types
        for k=1:numel(type_string_s)
            type_string=type_string_s{k};
            if isa(data, type_string)
                tp=k-1;
                break;
            end
        end
    end

    if isempty(tp)
        error('Data not understood');
    end

    vec_typ=repmat(tp,1,vec_num);


function f=get_ascii_printer(vec_typ,data)
% use vec_typ to create a format string for printf
% data is used in case we print floats; use the least possible number of
% floating point positions

    ni_defs=afni_ni_defs();

    n_vec_typ=numel(vec_typ);

    if vec_typ_and_data_is_all_numeric(vec_typ, data)
        pat1=get_single_element_ascii_print_format(vec_typ,data,ni_defs);
        pat=repmat([pat1 ' '],1,n_vec_typ);
        pat(end)=sprintf('\n');

        f=@(x)sprintf(pat,x');
    else
        if ~iscell(data)
            error('illegal data: espected cell, found %s', class(data));
        end

        n_col=numel(data);
        if n_vec_typ~=n_col
            error('vec_typ has %d elements, but data has %d elements',...
                        n_vec_typ, n+col);
        end

        pats=cell(1,n_col);
        for k=1:n_col
            pats{k}=get_single_element_ascii_print_format(vec_typ(k),...
                                                        data{k},ni_defs);
        end

        f=@(x)print_cell_ascii(pats,data);
    end

function tf=vec_typ_and_data_is_all_numeric(vec_typ, data)
    ni_defs=afni_ni_defs();
    tf=all(vec_typ(1)~=ni_defs.NI_STRING) && ~iscell(data);



function s=print_cell_ascii(pats, data_cell)
% mixed data types, or anything with a string
    string_mask=cellfun(@iscellstr,data_cell);
    n_col=numel(pats);

    n_row_s=cellfun(@numel,data_cell);

    n_row=n_row_s(1);
    if any(n_row~=n_row_s);
        error('data has different number of elements');
    end

    col_sep=' ';
    row_sep=sprintf('\n');

    s_cell=cell(1,n_row*n_col);
    pos=0;
    for row=1:n_row
        for col=1:n_col
            d_col=data_cell{col};

            if string_mask(col)
                d=niml_escape_string(d_col{row});
            else
                d=d_col(row);
            end

            pos=pos+1;

            if col==n_col
                sep=row_sep;
            else
                sep=col_sep;
            end

            s_cell{pos}=[sprintf(pats{col}, d) sep];

        end
    end

    s=[row_sep cat(2,s_cell{:})];



function s_escaped=niml_escape_string(s)
    ni_defs=afni_ni_defs();
    escape=ni_defs.escape;

    s_escaped=char(s);

    n=size(escape,1);
    for k=1:n
        s_escaped=strrep(s_escaped,escape{k,2},escape{k,1});
    end




function f=get_single_element_ascii_print_format(vec_typ,data,ni_defs)
    vec_str=strtrim(ni_defs.type_name(vec_typ+1,:));

    switch vec_str
        case 'String'
            f='"%s"';

        case {'byte','short','int'}
            % CHECKME don't know if this is ok for bytes and shorts...
            f='%d';

        otherwise
            precision=get_required_precision(data);
            f=sprintf('%%.%df',precision);
    end


