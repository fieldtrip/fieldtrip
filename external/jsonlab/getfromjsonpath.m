function obj=getfromjsonpath(root, jsonpath)
%
%    obj=getfromjsonpath(root, jsonpath)
%
%    Query and retrieve elements from matlab data structures using JSONPath
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        root: a matlab data structure like an array, cell, struct, etc
%        jsonpath: a string in the format of JSONPath, see loadjson help
%
%    output:
%        obj: if the specified element exist, obj returns the result
%
%    example:
%        getfromjsonpath(struct('a',[1,2,3]), '$.a.[1]')      % returns 2
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

obj=root;
[pat,paths]=regexp(jsonpath,'\.*([^\s\.]+)\.*','match','tokens');
if(~isempty(pat) && ~isempty(paths))
for i=1:length(paths)
    if(strcmp(paths{i}{1},'$'))
        continue;
    elseif(regexp(paths{i}{1},'$\d+'))
        obj=obj(str2double(paths{i}{1}(2:end))+1);
    elseif(regexp(paths{i}{1},'^\[\d+\]$'))
        if(iscell(obj))
            obj=obj{str2double(paths{i}{1}(2:end-1))+1};
        else
            obj=obj(str2double(paths{i}{1}(2:end-1))+1);
        end
    elseif(isstruct(obj))
        obj=obj.(encodevarname(paths{i}{1}));
    elseif(isa(obj,'containers.Map'))
        obj=obj(paths{i}{1});
    elseif(isa(obj,'table'))
        obj=obj(:,paths{i}{1});
    else
        error('json path segment (%d) "%s" can not be found in the input object\n',i,paths{i}{1});
    end
end
end