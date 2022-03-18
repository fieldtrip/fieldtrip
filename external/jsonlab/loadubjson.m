function varargout = loadubjson(varargin)
%
% data=loadubjson(fname,opt)
%    or
% data=loadubjson(fname,'param1',value1,'param2',value2,...)
%
% Parse a UBJSON file or string and store the output into a MATLAB variable
%
% author: Qianqian Fang (q.fang <at> neu.edu)
% initially created on 2019/06/08
%
% This function is an alias to loadbj
%
% input:
%      fname: input file name, if the file with such name exists, it will
%             be read, otherwise, this function will attempt to parse the
%             string in fname as a UBJSON stream
%      opt: a struct to store parsing options, opt can be replaced by 
%           a list of ('param',value) pairs - the param string is equivallent
%           to a field in opt. The supported options can be found by typing
%           "help loadbj".
%
% output:
%      data: a cell array, where {...} blocks are converted into cell arrays,
%           and [...] are converted to arrays
%
% examples:
%      obj=struct('string','value','array',[1 2 3]);
%      ubjdata=saveubjson('obj',obj);
%      dat=loadubjson(ubjdata)
%      dat=loadubjson(['examples' filesep 'example1.ubj'])
%      dat=loadubjson(['examples' filesep 'example1.ubj'],'SimplifyCell',0)
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

[varargout{1:nargout}]=loadbj(varargin{:},'endian','B');
