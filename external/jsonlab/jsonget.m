function json=jsonget(fname,mmap,varargin)
%
% json=jsonget(fname,mmap,'$.jsonpath1','$.jsonpath2',...)
%
% Fast reading of JSON data records using memory-map (mmap) returned by
% loadjson and JSONPath-like keys
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
% initially created on 2022/02/02
%
% input:
%      fname: a JSON/BJData/UBJSON string or stream, or a file name
%      mmap: memory-map returned by loadjson/loadbj of the same data
%            important: mmap must be produced from the same file/string,
%            otherwise calling this function may cause data corruption
%      '$.jsonpath1,2,3,...':  a series of strings in the form of JSONPath
%            as the key to each of the record to be retrieved; if no paths
%            are given, all items in mmap are retrieved
%
% output:
%      json: a cell array, made of elements {'$.jsonpath_i',json_string_i}
%
% examples:
%      str='[[1,2],"a",{"c":2}]{"k":"test"}';
%      [dat, mmap]=loadjson(str);
%      savejson('',dat,'filename','mydata.json','compact',1);
%      json=jsonget(str,mmap,'$.[0]','$.[2].c')
%      json=jsonget('mydata.json',mmap,'$.[0]','$.[2].c')
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

if(ischar(fname) || isa(fname,'string'))
        if(regexp(fname,'^\s*(?:\[.*\])|(?:\{.*\})\s*$','once'))
            inputstr=fname;
        elseif(~exist('memmapfile','file'))
            if(exist(fname,'file'))
               try
                   fid = fopen(fname,'rb');
               catch
               end
            end
        end
end

mmap=[mmap{:}];
keylist=mmap(1:2:end);

loc=1:length(keylist);
if(length(varargin)>=1)
    [tf,loc]=ismember(varargin,keylist);
    if(any(tf))
       keylist=keylist(loc);
    else
       keylist={};
    end
end

json={};

if(isstruct(fname) || iscell(fname) || isa(fname,'table') || isa(fname,'containers.Map'))
    for i=1:length(keylist)
        json{end+1}=getfromjsonpath(fname,keylist{i});
    end
    if(length(json)==1)
        json=json{1};
    end
    return;
end

for i=1:length(keylist)
    bmap=mmap{loc(i)*2};
    rec={'uint8',[1,bmap(2)],  'x'};
    if(exist('inputstr','var'))
        json{end+1}={keylist{i}, inputstr(bmap(1):bmap(1)+bmap(2)-1)};
    else
        if(exist('fid','var') && fid>=0)
            fseek(fid, bmap(1), -1);
            json{end+1}={keylist{i}, fread(fid,bmap(1),'uint8=>char')};
        else
            fmap=memmapfile(fname,'writable',false, 'offset',bmap(1),'format', rec);
            json{end+1}={keylist{i}, char(fmap.Data(1).x)};
        end
    end
end

if(exist('fid','var') && fid>0)
    fclose(fid);
end