function niml_cell = afni_niml_parse(s)
% Simple parsing routine for AFNI NIML datasets
%
% N=AFNI_NIML_PARSE(S) parses the niml string S and returns a parsed
% structure P.
%
% P can either be an niml element, or a group of niml elements. In the
% formder case, P contains a field .data; in the latter, a cell .nodes.
%
% This function is more or less the inverse of AFNI_NIML_PRINT.
%
% Many thanks to Ziad Saad, who wrote afni_niml_read and other routines
% that formed the basis of this function.
%
% Please note that this function is *VERY EXPERIMENTAL*
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>
% NNO May 2015 <nikolaas.oosterhof@unitn.it>

    len_s=numel(s);

    niml_cell=cell(1,1);

    pos=1;
    elem_counter=0;
    while pos<len_s
        [niml,pos]=niml_parse_helper(s,1);

        elem_counter=elem_counter+1;
        if numel(niml_cell)>elem_counter
            % allocate more space
            niml_cell{1+2*elem_counter,1}=[];
        end

        niml_cell{elem_counter}=niml;
    end

    niml_cell=niml_cell(1:elem_counter);


function [niml,pos]=niml_parse_helper(s, pos)
% pos is current character that has to be parsed
% i.e. pos=1 means the very beginning of the string

    pos=consume_whitespace_optionally(s, pos);

    if starts_with(s,pos,'<xml')
        % ignore anything in between xml characters
        pos=find_char(s,pos,'>');
        pos=consume_whitespace_optionally(s, pos);
    end

    header_start_pos=find_char(s, pos, '<');

    key_values_start_pos=find_whitespace(s, header_start_pos);

    header_name=char(s(header_start_pos:(key_values_start_pos-2)));
    header_end_pos=find_char(s, key_values_start_pos, '>');

    s_key_values=s(key_values_start_pos:(header_end_pos-2));

    niml=niml_parse_header_key_values(s_key_values);
    niml.name=header_name;

    assert(header_name(1)~='/');

    is_group=isfield(niml,'ni_form') && ...
                    strcmp(niml.ni_form,'ni_group');
    if is_group
        [niml,pos]=niml_parse_group(s, header_end_pos, niml);
    else
        [niml,pos]=niml_parse_data(s, header_end_pos, niml);
    end

    pos=consume_whitespace_optionally(s, pos);
    pos=niml_consume_endmarker(s, pos, niml);
    pos=consume_whitespace_optionally(s, pos);



function [niml,pos]=niml_parse_group(s, pos, niml)
    end_marker=['</' niml.name '>'];

    niml_nodes=cell(5,1);
    n_niml_nodes=0;
    while true
        pos=consume_whitespace_optionally(s, pos);

        if starts_with(s, pos, end_marker)
            niml.nodes=niml_nodes(1:n_niml_nodes);
            return
        end

        n_niml_nodes=n_niml_nodes+1;

        if n_niml_nodes>numel(niml_nodes)
            % double the space
            niml_nodes{1+2*n_niml_nodes}=[];
        end

        % use recursive call
        [niml_nodes{n_niml_nodes},pos]=niml_parse_helper(s, pos);
    end


function [niml,pos]=niml_parse_data(s, pos, niml)
% parse data in body
    niml.vec_typ = afni_nel_getvectype(niml.ni_type);
    niml.vec_len = str2num(niml.ni_dimen);
    niml.vec_num = length(niml.vec_typ);
    if afni_ni_is_numeric_type(niml.vec_typ)
        [niml,pos] = niml_parse_numeric_data(s, pos, niml);
    else
        [niml,pos]= niml_parse_mixed_data(s, pos, niml);
    end


function [niml,pos] = niml_parse_numeric_data(s, pos, niml)
% parse numeric data in body, supporting binary data
    n_columns=niml.vec_num;
    n_rows=niml.vec_len;

    is_binary=isfield(niml,'ni_form') && ...
            starts_with(niml.ni_form,1,'binary');

    if is_binary
        ni_def=afni_ni_defs();

        vec_typ=niml.vec_typ;
        if ~all(vec_typ(1)==vec_typ)
            error('Mixed data not supported for binary data');
        end

        ni_def_index=niml.vec_typ(1)+1; % base 1
        bytes_per_element=ni_def.size_bytes(ni_def_index);
        n_data_bytes=bytes_per_element*n_columns*n_rows;

        data_type=func2str(ni_def.type{ni_def_index});


        data_bytes=uint8(s(pos+(0:(n_data_bytes-1))));
        data=typecast(data_bytes,data_type);

        d_pos=n_data_bytes+1;

        if niml_binary_data_requires_byteswap(niml.ni_form)
            data=afni_swapbytes(data);
        end

    else
        % data in ASCII form
        data_string=char(s(pos:end));
        [data,unused,unused,d_pos]=sscanf(data_string,...
                                        '%f',n_columns*n_rows);
    end

    niml.data = reshape(double(data), niml.vec_num, niml.vec_len)';

    pos=pos+d_pos-1;


function tf=niml_binary_data_requires_byteswap(ni_form)
    [unused,unused,endian]=computer();
    switch lower(ni_form)
        case 'binary.lsbfirst'
            tf=endian~='L';

        case 'binary.msbfirst'
            tf=endian~='M';

        case 'binary'
            % native
            tf=false;

        otherwise
            error('unrecognized ni_form ''%s''', ni_form);
    end


function [niml, pos]=niml_parse_mixed_data(s, pos, niml)
% parse mixed data in ASCII format
    pos=consume_whitespace_optionally(s, pos);
    n_col=niml.vec_num;
    n_row=niml.vec_len;
    vec_typ=niml.vec_typ;

    col_is_numeric=true(1,n_col);

    row_data_cell=cell(1,n_col);

    for row=1:n_row
        is_first_row=row==1;

        for col=1:n_col

            [data_elem,pos]=niml_parse_string_data_element(s, pos, ...
                                                    vec_typ(col));

            row_data_cell{col}=data_elem;

            if is_first_row
                col_is_numeric(col)=isnumeric(data_elem);
            end
        end

        if is_first_row
            all_col_are_numeric=all(col_is_numeric);
            if all_col_are_numeric
                data=zeros(n_row,n_col);
            else
                data=cell(n_row,n_col);
            end
        end

        if all_col_are_numeric
            row_data=cat(2,row_data_cell{:});
        else
            row_data=row_data_cell;
        end

        data(row,:)=row_data;
    end

    if all(col_is_numeric)
        niml.data=data;
    else
        niml.data=cell(1,n_col);
        for col=1:n_col
            col_data=data(:,col);
            niml.data{col}=cat(1,col_data{:});
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p=niml_parse_header_key_values(s)
% parse a string of the form "K1=V1 K2=V2 ...
    expr='\s*(?<lhs>\w+)\s*=\s*"(?<rhs>[^"]+)"';
    hh=regexp(char(s),expr,'names');

    if numel(hh)==1
        % Octave regexp output
        lhs=hh.lhs;
        rhs=hh.rhs;
    else
        % Matlab regexp output
        lhs={hh(:).lhs};
        rhs={hh(:).rhs};
    end

    p = cell2struct(rhs, lhs,2); clear hh;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse helper functions

function vec_typ=afni_nel_getvectype(tt)
% gets the vector type for a data element
%
% TODO: this function needs refactoring
    vec_typ=zeros(1,1000);
    nn=0;
    while (~isempty(tt)),
        [ttt,tt] = strtok(tt, ',');
        %look for the N*type syntax
        N = 1;
        [tttt,ttt] = strtok(ttt,'*');
        Ntttt = str2double(tttt);
        if (~isnan(Ntttt)),  %have a number, get the type
            N = Ntttt;
            tttt = strtok(ttt,'*');
        end
        vec_typ(nn+1:1:nn+N) = afni_ni_rowtype_name_to_code(tttt);
        nn = nn + N;
    end

    vec_typ=vec_typ(1:nn)-1; % convert to base0, as in the niml.h.


function pos=niml_consume_endmarker(s, pos, niml)
    %end_marker=['</' niml.name '>'];
    pos=consume_string(s, pos, '</');
    pos=consume_whitespace_optionally(s, pos);
    pos=consume_string(s, pos, niml.name);
    pos=consume_string(s, pos, '>');

    %pos=consume_string(s, pos, end_marker);

function [data,pos]=niml_parse_string_data_element(s, pos, vec_typ)
    % parse data in ascii format

    n_s=numel(s);
    if pos >= n_s
        error('no data to parse');
    end

    % define delimeters
    quote=uint8('"''');
    data_delim_char=uint8(';');
    ni_def=afni_ni_defs();

    % process white space
    pos=consume_whitespace_optionally(s, pos);

    % if quote, process it and use quote as delimeter
    has_quote=any(s(pos)==quote);

    if has_quote
        pos=pos+1;
        delimeter_chars=quote;
    else
        delimeter_chars=[get_whitespace_characters() data_delim_char];
    end

    % start of data element in string representation
    start_pos=pos;

    % take all characters until delimeter is found
    while pos < n_s && ~any(s(pos)==delimeter_chars)
        pos=pos+1;
    end

    data_value=char(s(start_pos:(pos-1)));

    % if quote
    if has_quote
        assert(any(s(pos)==quote)); % enforced by code above
        pos=pos+1;
    end

    % process remaining whitespace
    pos=consume_whitespace_optionally(s, pos);

    string_def_index=strmatch('String',ni_def.type_name)-1;
    is_string=vec_typ==string_def_index;

    if is_string
        % convert data
        data={niml_unescape_string(data_value)};
    else
        % convert to numeric
        data=str2double(data_value);
        assert(numel(data)==1);
    end


function s_unescaped=niml_unescape_string(s)
    ni_def=afni_ni_defs();
    escape=ni_def.escape;

    s_unescaped=char(s);

    n=size(escape,1);
    for k=n:-1:1
        s_unescaped=strrep(s_unescaped,escape{k,1},escape{k,2});
    end



%%%%%%%%%%%%%%%%%%%%%%%%%
% general string matching
function pos=consume_string_optionally(s, pos, needle)
    % skips needle if s starts with it at pos
    if starts_with(s, pos, needle)
        pos=pos+numel(needle);
    end

function pos=consume_string(s, pos, needle)
    % skips needle if s starts with it at pos, raises an error if that is
    % not the case
    c_pos=consume_string_optionally(s, pos, needle);
    if c_pos==pos
        error('did not find ''%s'' at position %d: ''%s ...'')', ...
                    needle, pos, s(pos+(0:10)));
    end
    pos=c_pos;

function pos=consume_char_optionally(s, pos, c)
    % skip any characters in c from s starting at pos
    n=numel(s);
    c_bytes=uint8(c);
    while pos<=n && any(s(pos)==c_bytes)
        pos=pos+1;
    end


function next_pos=consume_whitespace_optionally(s, pos)
    % finds the next character position that is not whitespace in
    % string s starting at pos
    next_pos=consume_char_optionally(s, pos, get_whitespace_characters());


function tf=starts_with(s, pos, needle)
    % sees whether s starts with needle at pos
    n_s=numel(s);
    n_needle=numel(needle);

    tf=n_needle<=(n_s-pos+1) && ...
                    strcmp(char(s(pos+(0:(n_needle-1)))),needle);

function char_pos=find_char(s, pos, c, negate)
    % finds the next character position after the first occurence
    % of c in string s starting at pos. c can contain multiple characters,
    % in which case the first occurence of either is matched.
    % If negate==true, then the first non-occurence of c is returned.
    if nargin<4
        negate=false;
    end

    if negate
        offset=0;
    else
        offset=1;
    end

    c_bytes=uint8(c);

    n=numel(s);
    for start_pos=pos:n
        if any(s(start_pos)==c_bytes)==~negate
            char_pos=start_pos+offset;
            return;
        end
    end

    error('Character %s not found', c);


function next_pos=find_whitespace(s, pos)
    % finds the next character position after the first occurence
    % of whitespace in string s starting at pos.
    whitespace_characters=get_whitespace_characters();
    next_pos=find_char(s, pos, whitespace_characters, false);


function whitespace_characters=get_whitespace_characters()
    whitespace_characters=uint8(sprintf(' \t\r\n\f\v'));

