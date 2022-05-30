function varargout=loadjd(filename, varargin)
%
%    data=loadjd(inputfile)
%       or
%    [data, mmap]=loadjd(inputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Parse a hierarchical container data file, including JSON,
%    binary JSON (BJData/UBJSON/MessagePack) and HDF5, and output
%    the parsed data in a MATLAB/Octave data structure
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        inputfile: the input hierarchical container data file, supporting:
%                *.json,.jnii,.jdt,.jmsh,.jnirs: JSON/JData based data files, see https://neurojson.org/jdata/draft2
%                *.bjd,.bnii,.jdb,.bmsh,.bnirs: binary JData (BJData) files, see https://neurojson.org/bjdata/draft2
%                *.ubj: UBJSON-encoded files, see http://ubjson.org
%                *.msgpack: MessagePack-encoded files, see http://msgpack.org
%                *.h5,.hdf5,.snirf: HDF5 files, see https://www.hdfgroup.org/
%        options: (optional) for JSON/JData files, these are optional 'param',value pairs
%                supported by loadjson.m; for BJData/UBJSON/MessagePack files, these are
%                options supported by loadbj.m; for HDF5 files, these are options
%                supported by loadh5.m (part of EasyH5 toolbox, http://github.com/fangq/easyh5/)
%
%    output:
%        data: a structure (array) or cell (array) storing the hierarchical data
%              in the container data file
%        mmap: (optional) output storing the JSON/binary JSON memory-map table for fast
%              disk access. see help info for loadjson or loadbj for more details.
%
%    examples:
%        obj=struct('string','value','array',[1 2 3]);
%        savejd('obj', obj, 'datafile.json')
%        newobj=loadjd('datafile.json');
%
%    license:
%        BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
%
% -- this function is part of JSONLab toolbox (http://iso2mesh.sf.net/cgi-bin/index.cgi?jsonlab)
%

if(nargin<1)
    error('you must provide file name');
end

if(regexpi(filename,'\.json$|\.jnii$\.jdt$\.jdat$|\.jmsh$|\.jnirs$'))
    [varargout{1:nargout}]=loadjson(filename,varargin{:});
elseif(regexpi(filename,'\.bjd$|\.bnii$\.jdb$\.jbat$|\.bmsh$|\.bnirs$'))
    [varargout{1:nargout}]=loadbj(filename,varargin{:});
elseif(regexpi(filename,'\.ubj$'))
    [varargout{1:nargout}]=loadubjson(filename,varargin{:});
elseif(regexpi(filename,'\.msgpack$'))
    [varargout{1:nargout}]=loadmsgpack(filename,varargin{:});
elseif(regexpi(filename,'\.h5$|\.hdf5$|\.snirf$'))
    if(~exist('loadh5','file'))
        error('you must first install EasyH5 from http://github.com/fangq/easyh5/');
    end
    [varargout{1:nargout}]=loadh5(filename,varargin{:});
else
    error('file suffix must be one of .json,.jnii,.jdt,.jmsh,.jnirs,.bjd,.bnii,.jdb,.bmsh,.bnirs,.ubj,.msgpack,.h5,.hdf5,.snirf');
end
