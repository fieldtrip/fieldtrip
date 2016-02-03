function s = convert(tree,uid)
% XMLTREE/CONVERT Converter an XML tree in a structure
% 
% tree      - XMLTree object
% uid       - uid of the root of the subtree, if provided.
%             Default is root
% s         - converted structure
%__________________________________________________________________________
%
% Convert an XMLTree into a structure, when possible.
% When several identical tags are present, a cell array is used.
% The root tag is not saved in the structure.
% If provided, only the structure corresponding to the subtree defined
% by the uid UID is returned.
%__________________________________________________________________________
% Copyright (C) 2002-2015  http://www.artefact.tk/

% Guillaume Flandin
% $Id$

% Exemple:
% tree = '<a><b>field1</b><c>field2</c><b>field3</b></a>';
% toto = convert(xmltree(tree));
% <=> toto = struct('b',{{'field1', 'field3'}},'c','field2')

%error(nargchk(1,2,nargin));

% Get the root uid of the output structure
if nargin == 1
    % Get the root uid of the XML tree
    root_uid = root(tree);
else
    % Uid provided by user
    root_uid = uid;
end

% Initialize the output structure
% struct([]) should be used but this only works with Matlab 6
% So we create a field that we delete at the end
%s = struct(get(tree,root_uid,'name'),''); % struct([])
s = struct('deletedummy','');

%s = sub_convert(tree,s,root_uid,{get(tree,root_uid,'name')});
s = sub_convert(tree,s,root_uid,{});

s = rmfield(s,'deletedummy');

%==========================================================================
function s = sub_convert(tree,s,uid,arg)
    type = get(tree,uid,'type');
    switch type
        case 'element'
            child = children(tree,uid);
            l = {};
            ll = {};
            for i=1:length(child)
                if isfield(tree,child(i),'name')
                    ll = { ll{:}, get(tree,child(i),'name') };
                end
            end
            for i=1:length(child)
                if isfield(tree,child(i),'name')
                    name = get(tree,child(i),'name');
                    nboccur = sum(ismember(l,name));
                    nboccur2 = sum(ismember(ll,name));
                    l = { l{:}, name };
                    if nboccur || (nboccur2>1)
                        arg2 = { arg{:}, name, {nboccur+1} };
                    else
                        arg2 = { arg{:}, name};
                    end
                else
                    arg2 = arg;
                end
                s = sub_convert(tree,s,child(i),arg2);
            end
            if isempty(child)
                s = sub_setfield(s,arg{:},'');
            end
            %- saving attributes : does not work with <a t='q'>b</a>
            %- but ok with <a t='q'><c>b</c></a>
%             attrb = attributes(tree,'get',uid);     %-
%             if ~isempty(attrb)                      %-
%                 arg2 = {arg{:} 'attributes'};       %-
%                 s = sub_setfield(s,arg2{:},attrb);  %-
%             end                                     %-
        case 'chardata'
            s = sub_setfield(s,arg{:},get(tree,uid,'value'));
            %- convert strings into their numerical equivalent when possible
            %- e.g. string '3.14159' becomes double scalar 3.14159
%             v = get(tree,uid,'value');              %-
%             cv = str2num(v);                        %-
%             if isempty(cv)                          %-
%                 s = sub_setfield(s,arg{:},v);       %-
%             else                                    %-
%                 s = sub_setfield(s,arg{:},cv);      %-
%             end                                     %-
        case 'cdata'
            s = sub_setfield(s,arg{:},get(tree,uid,'value'));
        case 'pi'
            % Processing instructions are evaluated if possible
            app = get(tree,uid,'target');
            switch app
                case {'matlab',''}
                    s = sub_setfield(s,arg{:},eval(get(tree,uid,'value')));
                case 'unix'
                    s = sub_setfield(s,arg{:},unix(get(tree,uid,'value')));
                case 'dos'
                    s = sub_setfield(s,arg{:},dos(get(tree,uid,'value')));
                case 'system'
                    s = sub_setfield(s,arg{:},system(get(tree,uid,'value')));
                otherwise
                    try
                        s = sub_setfield(s,arg{:},feval(app,get(tree,uid,'value')));
                    catch
                        warning('[XMLTree] Unknown target application');
                    end
            end
        case 'comment'
            % Comments are forgotten
        otherwise
            warning(sprintf('Type %s unknown : not saved',get(tree,uid,'type')));
    end
    
%==========================================================================
function s = sub_setfield(s,varargin)
% Same as setfield but using '{}' rather than '()'
%if (isempty(varargin) | length(varargin) < 2)
%    error('Not enough input arguments.');
%end

subs = varargin(1:end-1);
for i = 1:length(varargin)-1
    if (isa(varargin{i}, 'cell'))
        types{i} = '{}';
    elseif ischar(varargin{i})
        types{i} = '.';
        subs{i} = varargin{i}; %strrep(varargin{i},' ',''); % deblank field name
    else
        error('Inputs must be either cell arrays or strings.');
    end
end

% Perform assignment
try
   s = builtin('subsasgn', s, struct('type',types,'subs',subs), varargin{end});
catch
   error(lasterr)
end
