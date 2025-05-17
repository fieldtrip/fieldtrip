function data = regrouph5(root, varargin)
%
%    data=regrouph5(root)
%       or
%    data=regrouph5(root,type)
%    data=regrouph5(root,{'nameA','nameB',...})
%
%    Processing a loadh5 restored data and merge "indexed datasets", whose
%    names start with an ASCII string followed by a contiguous integer
%    sequence number starting from 1, into a cell array. For example,
%    datasets {data.a1, data.a2, data.a3} will be merged into a cell/struct
%    array data.a with 3 elements.
%
%    A single subfield .name1 will be renamed as .name. Items with
%    non-contigous numbering will not be grouped. If .name and
%    .name1/.name2 co-exist in the input struct, no grouping will be done.
%
%    The grouped subfield will appear at the position of the first
%    pre-grouped item in the original input structure.
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        root: the raw input HDF5 data structure (loaded from loadh5.m)
%        type: if type is set as a cell array of strings, it restrict the
%              grouping only to the subset of field names in this list;
%              if type is a string as 'snirf', it is the same as setting
%              type as {'aux','data','nirs','stim','measurementList'}.
%
%    output:
%        data: a reorganized matlab structure.
%
%    example:
%        a=struct('a1',rand(5),'a2','string','a3',true,'d',2+3i,'e',{'test',[],1:5});
%        regrouph5(a)
%        saveh5(a,'test.h5');
%        rawdata=loadh5('test.h5')
%        data=regrouph5(rawdata)
%
%    this file is part of EasyH5 Toolbox: https://github.com/NeuroJSON/easyh5
%
%    License: GPLv3 or 3-clause BSD license, see https://github.com/NeuroJSON/easyh5 for details
%

if (nargin < 1)
    help regrouph5;
    return
end

dict = {};
if (~isempty(varargin))
    if (ischar(varargin{1}) && strcmpi(varargin{1}, 'snirf'))
        dict = {'aux', 'data', 'nirs', 'stim', 'measurementList'};
    elseif (iscell(varargin{1}))
        dict = varargin{1};
    end
end

data = struct;
if (isstruct(root))
    data = repmat(struct, size(root));
    names = fieldnames(root);
    newnames = struct();
    firstpos = struct();

    for i = 1:length(names)
        item = regexp(names{i}, '^(.*\D)(\d+)$', 'tokens');
        if (~isempty(item) && str2double(item{1}{2}) ~= 0 && ~isfield(root, item{1}{1}))
            if (~isfield(newnames, item{1}{1}))
                newnames.(item{1}{1}) = str2double(item{1}{2});
            else
                newnames.(item{1}{1}) = [newnames.(item{1}{1}), str2double(item{1}{2})];
            end
            if (~isfield(firstpos, item{1}{1}))
                firstpos.(item{1}{1}) = length(fieldnames(data(1)));
            end
        else
            for j = 1:length(root)
                if (isstruct(root(j).(names{i})))
                    data(j).(names{i}) = regrouph5(root(j).(names{i}));
                else
                    data(j).(names{i}) = root(j).(names{i});
                end
            end
        end
    end

    names = fieldnames(newnames);
    if (~isempty(dict))
        names = intersect(names, dict);
    end

    for i = length(names):-1:1
        len = length(newnames.(names{i}));
        idx = newnames.(names{i});
        if ((min(idx) ~= 1 || max(idx) ~= len) && len ~= 1)
            for j = 1:len
                dataname = sprintf('%s%d', names{i}, idx(j));
                for k = 1:length(root)
                    if (isstruct(root(k).(dataname)))
                        data(k).(dataname) = regrouph5(root(k).(dataname));
                    else
                        data(k).(dataname) = root(k).(dataname);
                    end
                end
            end
            pos = firstpos.(names{i});
            len = length(fieldnames(data));
            data = orderfields(data, [1:pos, len, pos + 1:len - 1]);
            continue
        end
        for j = 1:length(data)
            data(j).(names{i}) = cell(1, len);
        end
        idx = sort(idx);
        for j = 1:len
            for k = 1:length(root)
                obj = root(k).(sprintf('%s%d', names{i}, idx(j)));
                if (isstruct(obj))
                    data(k).(names{i}){j} = regrouph5(obj);
                else
                    data(k).(names{i}){j} = obj;
                end
            end
        end
        pos = firstpos.(names{i});
        len = length(fieldnames(data));
        data = orderfields(data, [1:pos, len, pos + 1:len - 1]);
        try
            data.(names{i}) = cell2mat(data.(names{i}));
        catch
        end
    end
end
