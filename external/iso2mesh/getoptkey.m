function val=getoptkey(key,default,varargin)
%
% val=getoptkey(key,default,opt)
%    or
% val=getoptkey(key,default,'key1',val1,'key2',val2, ...)
%
% query the value of a key from a structure or a list of key/value pairs
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%   key: a string name for the target struct field name
%   default: the default value of the key is not found
%   opt: a struct object; the field names will be searched to match the 
%        key input, opt can be a list of 'keyname'/value pairs
%
% output:
%   val: val=opt.key if found, otherwise val=default
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

val=default;
if(nargin<=2) return; end
opt=varargin2struct(varargin{:});
if(isstruct(opt) && isfield(opt,key))
    val=getfield(opt,key);
end

