function varargout = blosc2encode(varargin)
%
% output = blosc2encode(input)
%    or
% [output, info] = blosc2encode(input)
%
% Compress a string or a numerical array using LZ4-compression
%
% This function depends on the ZMat toolbox (http://github.com/NeuroJSON/zmat)
%
% authors:Qianqian Fang (q.fang <at> neu.edu)
%
% input:
%      input: the original data, can be a string, a numerical vector or array
%
% output:
%      output: the compressed byte stream stored in a uint8 vector
%      info: (optional) a struct storing the metadata of the input, see "help zmat" for details
%
% examples:
%      [bytes, info]=blosc2encode(eye(10));
%      orig=blosc2decode(bytes,info);
%
% license:
%     BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

if (nargin == 0)
    error('you must provide at least 1 input');
end

if (exist('zmat', 'file') == 2 || exist('zmat', 'file') == 3)
    if (nargin > 1 && ischar(varargin{2}))
        [varargout{1:nargout}] = zmat(varargin{1}, 1, varargin{2:end});
    else
        [varargout{1:nargout}] = zmat(varargin{1}, 1, 'blosc2blosclz', varargin{2:end});
    end
    return
else
    error('you must install ZMat toolbox to use this feature: http://github.com/NeuroJSON/zmat');
end
