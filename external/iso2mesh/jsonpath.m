function obj = jsonpath(root, jpath, varargin)
%
%    obj=jsonpath(root, jpath)
%
%    Query and retrieve elements from matlab data structures using JSONPath
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        root: a matlab data structure like an array, cell, struct, etc
%        jpath: a string in the format of JSONPath, see loadjson help
%
%    output:
%        obj: if the specified element exist, obj returns the result
%
%    example:
%        jsonpath(struct('a',[1,2,3]), '$.a[1]')      % returns 2
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

obj = root;
jpath = regexprep(jpath, '([^.\]])(\[[-0-9:\*]+\])', '$1.$2');
jpath = regexprep(jpath, '\[[''"]*([^\]''"]+)[''"]*\]', '.[$1]');
jpath = regexprep(jpath, '\\\.', '_0x2E_');
while (regexp(jpath, '(\[[''"]*[^\]''"]+)\.(?=[^\]''"]+[''"]*\])'))
    jpath = regexprep(jpath, '(\[[''"]*[^\]''"]+)\.(?=[^\]''"]+[''"]*\])', '$1_0x2E_');
end
[pat, paths] = regexp(jpath, '(\.{0,2}[^\.]+)', 'match', 'tokens');
if (~isempty(pat) && ~isempty(paths))
    if (strcmp(paths{1}, '$'))
        paths(1) = [];
    end
    for i = 1:length(paths)
        [obj, isfound] = getonelevel(obj, paths, i, varargin{:});
        if (~isfound)
            return
        end
    end
end

%% scan function

function [obj, isfound] = getonelevel(input, paths, pathid, varargin)

pathname = paths{pathid};
if (iscell(pathname))
    pathname = pathname{1};
end

deepscan = ~isempty(regexp(pathname, '^\.\.', 'once'));

origpath = pathname;
pathname = regexprep(pathname, '^\.+', '');

if (strcmp(pathname, '$'))
    obj = input;
elseif (regexp(pathname, '$\d+'))
    obj = input(str2double(pathname(2:end)) + 1);
elseif (~isempty(regexp(pathname, '^\[[-0-9\*:]+\]$', 'once')) || iscell(input))
    arraystr = pathname(2:end - 1);
    if (find(arraystr == ':'))
        arrayrange = regexp(arraystr, '(?<start>-*\d*):(?<end>-*\d*)', 'names');
        if (~isempty(arrayrange.start))
            arrayrange.start = str2double(arrayrange.start);
            arrayrange.start = (arrayrange.start < 0) * length(input) + arrayrange.start + 1;
        else
            arrayrange.start = 1;
        end
        if (~isempty(arrayrange.end))
            arrayrange.end = str2double(arrayrange.end);
            arrayrange.end = (arrayrange.end < 0) * length(input) + arrayrange.end + 1;
        else
            arrayrange.end = length(input);
        end
    elseif (regexp(arraystr, '^[-0-9:]+', 'once'))
        firstidx = str2double(arraystr);
        if (firstidx < 0)
            firstidx = length(input) + firstidx + 1;
        else
            firstidx = firstidx + 1;
        end
        arrayrange.start = firstidx;
        arrayrange.end = firstidx;
    elseif (~isempty(regexp(arraystr, '^\*', 'once')))
        % do nothing
    end
    if (exist('arrayrange', 'var'))
        obj = {input(arrayrange.start:arrayrange.end)};
    else
        arrayrange = struct('start', 1, 'end', length(input));
    end
    if (~exist('obj', 'var') && iscell(input))
        input = {input{arrayrange.start:arrayrange.end}};
        if (deepscan)
            searchkey = ['..' pathname];
        else
            searchkey = origpath;
        end
        for idx = 1:length(input)
            [val, isfound] = getonelevel(input{idx}, [paths{1:pathid - 1} {searchkey}], pathid);
            if (isfound)
                if (~exist('newobj', 'var'))
                    newobj = {};
                end
                newobj = [newobj(:)', val(:)'];
            end
        end
        if (exist('newobj', 'var'))
            obj = newobj;
        end
        if (exist('obj', 'var') && iscell(obj) && length(obj) == 1)
            obj = obj{1};
        end
    else
        obj = input(arrayrange.start:arrayrange.end);
    end
elseif (isstruct(input) || isa(input, 'containers.Map') || isa(input, 'table'))
    pathname = regexprep(pathname, '^\[(.*)\]$', '$1');
    stpath = encodevarname(pathname);
    if (isstruct(input))
        if (isfield(input, stpath))
            obj = {input.(stpath)};
        end
    elseif (isa(input, 'table'))
        if (any(ismember(input.Properties.VariableNames, stpath)))
            obj = {input.(stpath)};
        end
    else
        if (isKey(input, pathname))
            obj = {input(pathname)};
        end
    end
    if (~exist('obj', 'var') || deepscan)
        if (isa(input, 'containers.Map'))
            items = keys(input);
        else
            items = fieldnames(input);
        end
        for idx = 1:length(items)
            if (isa(input, 'containers.Map'))
                [val, isfound] = getonelevel(input(items{idx}), [paths{1:pathid - 1} {['..' pathname]}], pathid);
            elseif (isa(input, 'table'))
                [val, isfound] = getonelevel({input.(items{idx})}, [paths{1:pathid - 1} {['..' pathname]}], pathid);
            elseif (isstruct(input) && length(input) > 1) % struct array
                [val, isfound] = getonelevel({input.(items{idx})}, [paths{1:pathid - 1} {['..' pathname]}], pathid);
            else
                [val, isfound] = getonelevel(input.(items{idx}), [paths{1:pathid - 1} {['..' pathname]}], pathid);
            end
            if (isfound)
                if (~exist('obj', 'var'))
                    obj = {};
                end
                if (iscell(val))
                    obj = [obj(:)', val(:)'];
                else
                    obj = [obj(:)', {val}];
                end
            end
        end
        if (exist('obj', 'var') && length(obj) == 1)
            obj = obj{1};
        end
    end
    if (exist('obj', 'var') && iscell(obj) && length(obj) == 1)
        obj = obj{1};
    end
elseif (~deepscan)
    error('json path segment "%s" can not be found in the input object\n', pathname);
end

if (~exist('obj', 'var'))
    isfound = false;
    obj = [];
elseif (nargout > 1)
    isfound = true;
end
