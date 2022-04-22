function ubj=saveubjson(rootname,obj,varargin)
%
% ubj=saveubjson(obj)
%    or
% ubj=saveubjson(rootname,obj,filename)
% ubj=saveubjson(rootname,obj,opt)
% ubj=saveubjson(rootname,obj,'param1',value1,'param2',value2,...)
%
% Convert a MATLAB object  (cell, struct, array, table, map, graphs ...) 
% into a Universal Binary JSON (UBJSON, Draft-12) or a MessagePack binary stream
%
% author: Qianqian Fang (q.fang <at> neu.edu)
% initially created on 2013/08/17
%
% Format specifications:
%    Binary JData (BJData):https://github.com/NeuroJSON/bjdata
%    UBJSON:               https://github.com/ubjson/universal-binary-json
%    MessagePack:          https://github.com/msgpack/msgpack
%
% This function is the same as calling "savebj(..,'ubjson',1,'endian','B')"
% By default this function creates UBJSON-compliant output without the
% newly added uint16(u), uint32(m), uint64(M) and half-precision float (h)
% data types and use Big-Endian for all numerical values as in UBJSON
% Draft-12.
%
% This function by default still enables an optimized ND-array format for efficient  
% array storage. To ensure the output compatible to UBJSON Draft-12, one should use
% "saveubjson(...,'NestArray',1)" or "savebj(...,'ubjson',1,'NestArray',1)"
%
% input:
%      rootname: the name of the root-object, when set to '', the root name
%           is ignored, however, when opt.ForceRootName is set to 1 (see below),
%           the MATLAB variable name will be used as the root name.
%      obj: a MATLAB object (array, cell, cell array, struct, struct array,
%           class instance)
%      filename: a string for the file name to save the output UBJSON data
%      opt: a struct for additional options, ignore to use default values.
%           opt can have the following fields (first in [.|.] is the default)
%
%           opt can be replaced by a list of ('param',value) pairs. The param 
%           string is equivallent to a field in opt and is case sensitive.
%
%           Please type "help savebj" for details for all supported options.
%
% output:
%      ubj: a binary string in the UBJSON format (see http://ubjson.org)
%
% examples:
%      jsonmesh=struct('MeshVertex3',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 
%               'MeshTet4',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],...
%               'MeshTri3',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;...
%                          2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],...
%               'MeshCreator','FangQ','MeshTitle','T6 Cube',...
%               'SpecialData',[nan, inf, -inf]);
%      saveubjson(jsonmesh)
%      saveubjson('',jsonmesh,'meshdata.ubj')
%      saveubjson('mesh1',jsonmesh,'FileName','meshdata.msgpk','MessagePack',1)
%      saveubjson('',jsonmesh,'KeepType',1)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

if(nargin==1)
    ubj=savebj('',rootname,'ubjson',1,'endian','B');
elseif(length(varargin)==1 && ischar(varargin{1}))
    ubj=savebj(rootname,obj,'FileName',varargin{1},'ubjson',1,'endian','B');
else
    ubj=savebj(rootname,obj,varargin{:},'ubjson',1,'endian','B');
end
