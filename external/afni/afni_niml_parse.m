function niml = afni_niml_parse(s)
% Simple parsing routine for AFNI NIML datasets
%
% N=AFNI_NIML_PARSE(S) parses the niml string S and returns a parsed
% structure P.
%
% P can either be an niml element, or a group of niml elements. In the
% formder case, P contains a field .data; in the latter, a cell .nodes.
%
% If S contains multiple niml datasets that are not in a group, then N will
% be a cell with the parsed datasets.
%
% This function is more or less the inverse of AFNI_NIML_PRINT. 
%
% Many thanks to Ziad Saad, who wrote afni_niml_read and other routines
% that formed the basis of this function.
%
% Please note that this function is *VERY EXPERIMENTAL*
%
% NNO Dec 2009 <n.oosterhof@bangor.ac.uk>

% if the input is a string (i.e. not called recursively), parse header and
% body and put them in a struct s.
if ischar(s)
    s=afni_nel_parseheaderbody(s);
end

% make sure we only have the fieldnames we expect
if ~isempty(setxor(fieldnames(s),{'headername','headertext','body'}))
    error('Illegal struct s, expected s.headername, s.headertext, s.body');
end

% if more than 1 element, parse each of them separately and put them in a
% cell.
scount=numel(s);
if scount>1
    niml=cell(scount,1);
    for k=1:scount
        niml{k}=afni_niml_parse(s(k));
    end
    return
end

% parse the header part
niml=afni_nel_parsekeyvalue(s.headertext);
niml.name=s.headername;

if isfield(niml,'ni_form') && strcmp(niml.ni_form,'ni_group')
   % this is a group of elements. parse each of the elements in the group 
   % and put the results in a field .nodes 
   
   nds=afni_niml_parse(s.body);
   
   % for consistency, ensure that .nodes is always a cell
   if iscell(nds)
       niml.nodes=nds;
   else
       niml.nodes{1}=nds;
   end
else
    % this is a single element
    
    % set a few fields
    niml.vec_typ = afni_nel_getvectype(niml.ni_type);
    niml.vec_len = str2num(niml.ni_dimen);
    niml.vec_num = length(niml.vec_typ);

    % parse only 
    if (~afni_ni_is_numeric_type(niml.vec_typ)),
      %fprintf(2,'Data not all numeric, will not parse it');
 
        niml.data = afni_nel_parse_nonnumeric(niml, s.body);
    else
        niml.data = afni_nel_parse_data(niml, s.body);
    end
end   
   
function p=afni_nel_parsekeyvalue(s)
% parses a string of the form "K1=V1 K2=V2 ...
    expr='\s*(?<lhs>\w+)\s*=\s*"(?<rhs>[^"]+)"';
    hh=regexp(s,expr,'names');
    p = cell2struct({hh(:).rhs}, {hh(:).lhs},2); clear hh;

function p=afni_nel_parseheaderbody(s)
% parses a header and body
% in the form <HEADERNAME HEADERTEXT>BODY</HEADERNAME>
    expr = '<(?<headername>\w+)(?<headertext>.*?)>(?<body>.*?)</\1>';
    p = regexp(s, expr,'names');
   
   
function vec_typ=afni_nel_getvectype(tt)   
% gets the vector type for a data element
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
                             %
                             % This is a point of concern as the
                             % afni_ni_rowtype_name_to_code function
                             % seems to prefer base1
                             %
                             % (if only matlab used base0 indexing...)
   

function p = afni_nel_parse_data(nel, data)
    d=sscanf(data,'%f');
    p = reshape(d, nel.vec_num, nel.vec_len)'; 
    
function p =afni_nel_parse_nonnumeric(nel,data)
    if strcmp(nel.ni_type,'String')
        p=strtrim(data);
        if strcmp(p([1 end]),'""')
            p=p(2:(end-1)); % remove surrounding quotes
        end
    else
        p=data; %do nothing
    end
    
        
   
