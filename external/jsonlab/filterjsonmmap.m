function mmap=filterjsonmmap(mmap, patterns, isinclude)
%
% mmap=filterjsonmmap(mmap, patterns, isinclude)
%
% filter JSON mmap keys based on inclusive or exclusive string patterns
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
% initially created on 2022/02/13
%
% input:
%      mmap: memory-map returned by loadjson/loadbj of the same data
%            important: mmap must be produced from the same file/string,
%            otherwise calling this function may cause data corruption
%      patterns: a string or a cell array of strings, each string will 
%            be tested to match the JSONPath keys in mmap
%      isinclude: 1 (default) to include all mmap entries that match at
%            least one of the patterns, and 0 - exclude those that match
%
% output:
%      mmap: a filtered JSON mmap
%
% examples:
%      str='{"arr":[[1,2],"a",{"c":2}],"obj":{"k":"test"}}';
%      [dat, mmap]=loadjson(str);
%      savejson('',mmap)
%      newmmap=filterjsonmmap(mmap,{'arr.[1]', 'obj.k'});
%      savejson('',newmmap)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%
    if(nargin<3)
        isinclude=1;
    end
    if(nargin>1 && ~isempty(patterns))
        keylist=[mmap{:}];
        keylist=keylist(1:2:end);
        if(~iscell(patterns))
            patterns={patterns};
        end
        mask=zeros(1,length(keylist));
        for i=1:length(patterns)
            mask=mask+cellfun(@length, strfind(keylist,patterns{i}));
        end
        if(isinclude)
            mmap=mmap(mask>0);
        else
            mmap(mask>0)=[];
        end
    end
end