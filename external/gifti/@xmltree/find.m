function list = find(varargin)
% XMLTREE/FIND Find elements in a tree with specified characteristics
% FORMAT list = find(varargin)
%
% tree  - XMLTree object
% xpath - string path with specific grammar (XPath)
% uid   - lists of root uid's
% parameter/value - pair of pattern
% list  - list of uid's of matched elements
%
%      list = find(tree,xpath)
%      list = find(tree,parameter,value[,parameter,value])
%      list = find(tree,uid,parameter,value[,parameter,value])
%
% Grammar for addressing parts of an XML document: 
%        XML Path Language XPath (http://www.w3.org/TR/xpath)
% Example: /element1//element2[1]/element3[5]/element4
%__________________________________________________________________________
%
% Find elements in an XML tree with specified characteristics or given
% a path (using a subset of XPath language).
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id$

% TODO:
%   - clean up subroutines
%   - find should only be documented using XPath (other use is internal)
%   - handle '*', 'last()' in [] and '@'
%   - if a key is invalid, should rather return [] than error ?

if nargin==0
    error('[XMLTree] A tree must be provided');
elseif nargin==1
    list = 1:length(tree.tree);
    return
elseif mod(nargin,2)
    list = sub_find_subtree1(varargin{1}.tree,root(varargin{1}),varargin{2:end});
elseif isa(varargin{2},'double') && ...
       ndims(varargin{2}) == 2 && ...
       min(size(varargin{2})) == 1
    list = unique(sub_find_subtree1(varargin{1}.tree,varargin{2:end}));
elseif nargin==2 && ischar(varargin{2})
    list = sub_pathfinder(varargin{:});
else
   error('[XMLTree] Arguments must be parameter/value pairs.');
end

%==========================================================================
function list = sub_find_subtree1(varargin)
    list = [];
    for i=1:length(varargin{2})
        res = sub_find_subtree2(varargin{1},...
                varargin{2}(i),varargin{3:end});
        list = [list res];
    end

%==========================================================================
function list = sub_find_subtree2(varargin)
    uid = varargin{2};
    list = [];
    if sub_comp_element(varargin{1}{uid},varargin{3:end})
        list = [list varargin{1}{uid}.uid];
    end
    if isfield(varargin{1}{uid},'contents')
        list = [list sub_find_subtree1(varargin{1},...
                varargin{1}{uid}.contents,varargin{3:end})];
    end

%==========================================================================
function match = sub_comp_element(varargin)
match = 0;
try
    % v = getfield(varargin{1}, varargin{2}); % slow...
    for i=1:floor(nargin/2)
        v = subsref(varargin{1}, struct('type','.','subs',varargin{i+1}));     
        if strcmp(v,varargin{i+2})
            match = 1;
        end
    end
catch
end

% More powerful but much slower
%match = 0;
%for i=1:length(floor(nargin/2)) % bug: remove length !!!
%   if isfield(varargin{1},varargin{i+1})
%     if ischar(getfield(varargin{1},varargin{i+1})) & ischar(varargin{i+2})
%      if strcmp(getfield(varargin{1},varargin{i+1}),char(varargin{i+2}))
%           match = 1;
%        end
%     elseif isa(getfield(varargin{1},varargin{i+1}),'double') & ...
%           isa(varargin{i+2},'double')
%        if getfield(varargin{1},varargin{i+1}) == varargin{i+2}
%           match = 1;
%       end
%    else 
%       warning('Cannot compare different objects');
%     end
%   end
%end

%==========================================================================
function list = sub_pathfinder(tree,pth)
    %- Search for the delimiter '/' in the path
    i = strfind(pth,'/');
    %- Begin search by root
    list = root(tree);
    %- Walk through the tree
    j = 1;
    while j <= length(i)
        %- Look for recursion '//'
        if j<length(i) && i(j+1)==i(j)+1
            recursive = 1;
            j = j + 1;
        else
            recursive = 0;
        end
        %- Catch the current tag 'element[x]'
        if j ~= length(i)
            element = pth(i(j)+1:i(j+1)-1);
        else
            element = pth(i(j)+1:end);
        end
        %- Search for [] brackets
        k = xml_findstr(element,'[',1,1);
        %- If brackets are present in current element
        if ~isempty(k)
            l   = xml_findstr(element,']',1,1);
            val = str2num(element(k+1:l-1));
            element = element(1:k-1);
        end
        %- Use recursivity
        if recursive
            list = find(tree,list,'name',element);
        %- Just look at children
        else
            if i(j)==1 % if '/root/...' (list = root(tree) in that case)
                if sub_comp_element(tree.tree{list},'name',element)
                    % list = 1; % list still contains root(tree)
                else
                    list = [];
                    return;
                end
            else
                list = sub_findchild(tree,list,element);
            end
        end
        % If an element is specified using a key
        if ~isempty(k)
            if val < 1 || val > length(list)+1
                error('[XMLTree] Bad key in the path.');
            elseif val == length(list)+1
                list = [];
                return;
            end
            list = list(val);
        end
        if isempty(list), return; end
        j = j + 1;
    end
    
%==========================================================================
function list = sub_findchild(tree,listt,elmt)
    list = [];
    for a=1:length(listt)
        for b=1:length(tree.tree{listt(a)}.contents)
            if sub_comp_element(tree.tree{tree.tree{listt(a)}.contents(b)},'name',elmt)
                list = [list tree.tree{tree.tree{listt(a)}.contents(b)}.uid];
            end
        end
    end
