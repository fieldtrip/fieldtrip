function json=jsonset(fname,mmap,varargin)
%
% json=jsonset(fname,mmap,'$.jsonpath1',newval1,'$.jsonpath2','newval2',...)
%
% Fast writing of JSON data records to stream or disk using memory-map 
% (mmap) returned by loadjson/loadbj and JSONPath-like keys
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
%            as the key to each of the record to be written
%
% output:
%      json: the modified JSON string or, in the case fname is a filename,
%            the cell string made of jsonpaths that are successfully
%            written
%
% examples:
%      % create test data
%       d.arr={[1,2],'a',struct('c',2)}; d.obj=struct('k','test') 
%      % convert to json string
%       str=savejson('',d,'compact',1)
%      % parse and return mmap
%       [dat, mmap]=loadjson(str);
%      % display mmap entries
%       savejson('',mmap)
%      % replace value using mmap
%       json=jsonset(str,mmap,'$.arr.[2].c','5')
%      % save same json string to file (must set savebinary 1)
%       savejson('',d,'filename','file.json','compact',1,'savebinary',1);
%      % fast write to file
%       json=jsonset('file.json',mmap,'$.arr.[2].c','5')
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

if(regexp(fname,'^\s*(?:\[.*\])|(?:\{.*\})\s*$','once'))
    inputstr=fname;
else
    if(~exist('memmapfile','file'))
        fid=fopen(fname,'r+b');
    end
end

mmap=[mmap{:}];
keylist=mmap(1:2:end);

opt=struct;
for i=1:2:length(varargin)
    if(isempty(regexp(varargin{i},'^\$', 'once')))
        opt.(encodevarname(varargin{i}))=varargin{i+1};
    end
end

json={};
for i=1:2:length(varargin)
    if(regexp(varargin{i},'^\$'))
        [tf,loc]=ismember(varargin{i},keylist);
        if(tf)
            bmap=mmap{loc*2};
            if(ischar(varargin{i+1}))
                val=varargin{i+1};
            else
                val=savejson('',varargin{i+1},'compact',1);
            end
            if(length(val)<=bmap(2))
                val=[val repmat(' ',[1,bmap(2)-length(val)])];
                if(exist('inputstr','var'))
                    inputstr(bmap(1):bmap(1)+bmap(2)-1)=val;
                else
                    if(exist('memmapfile','file'))
                        rec={'uint8', [1 bmap(2)],  'x'};
                        fmap=memmapfile(fname,'writable',true,'offset', bmap(1)-1, 'format', rec, 'repeat',1);
                        fmap.Data.x=uint8(val);
                    else
                        fseek(fid,bmap(1)-1,'bof');
                        fwrite(fid,val);
                    end
                    json{end+1}={varargin{i},val};
                end
            end
        end
    end
end

if(exist('fid','var') && fid>=0)
    fclose(fid);
end

if(exist('inputstr','var'))
    json=inputstr;
end
