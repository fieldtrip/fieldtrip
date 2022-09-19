function varargout=savejd(varargin)
%
%    savejd(rootname, obj, outputfile)
%       or
%    buffer=savejd(obj)
%    buffer=savejd(rootname, obj)
%    savejd(rootname, obj, opt)
%    savejd(rootname, obj, 'filename', outputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Save a complex MATLAB/Octave data structure to a hierarchical
%    container data file including JSON, binary JSON
%    (BJData/UBJSON/MessagePack) and HDF5
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        rootname: the name of the root-object, when set to '', the root name
%           is ignored, however, when opt.ForceRootName is set to 1 (see help savejson),
%           the MATLAB variable name will be used as the root name.
%        obj: a MATLAB object (array, cell, cell array, struct, struct array,
%             class instance).
%        outputfile: the output file name to the hierarchical container file
%                *.json,.jnii,.jdt,.jmsh,.jnirs: JSON/JData based data files, see https://neurojson.org/jdata/draft2
%                *.bjd,.bnii,.jdb,.bmsh,.bnirs: binary JData (BJData) files, see https://neurojson.org/bjdata/draft2
%                *.ubj: UBJSON-encoded files, see http://ubjson.org
%                *.msgpack: MessagePack-encoded files, see http://msgpack.org
%                *.h5,.hdf5,.snirf: HDF5 files, see https://www.hdfgroup.org/
%        opt: (optional) for JSON/JData files, these are optional 'param',value pairs
%                supported by loadjson.m; for BJData/UBJSON/MessagePack files, these are
%                options supported by loadbj.m; for HDF5 files, these are options
%                supported by loadh5.m (part of EasyH5 toolbox, http://github.com/fangq/easyh5/)
%
%    output:
%        data: a structure (array) or cell (array) storing the hierarchical data
%              in the container data file
%        mmap: optional output storing the JSON/binary JSON memory-map table for fast
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

opt=struct;

if(nargin>2)
    if(nargin==3 && ischar(varargin{3}))
        filename=varargin{3};
    else
        opt=varargin2struct(varargin{3:end});
        filename=jsonopt('filename','.json',opt);
    end
end
if(~exist('filename','var'))
    filename='.json';
end
if(regexpi(filename,'\.json$|\.jnii$\.jdt$\.jdat$|\.jmsh$|\.jnirs$'))
    [varargout{1:nargout}]=savejson(varargin{:});
elseif(regexpi(filename,'\.bjd$|\.bnii$\.jdb$\.jbat$|\.bmsh$|\.bnirs$'))
    [varargout{1:nargout}]=savebj(varargin{:});
elseif(regexpi(filename,'\.ubj$'))
    [varargout{1:nargout}]=saveubjson(varargin{:});
elseif(regexpi(filename,'\.msgpack$'))
    [varargout{1:nargout}]=savemsgpack(varargin{:});
elseif(regexpi(filename,'\.h5$|\.hdf5$|\.snirf$'))
    if(~exist('saveh5','file'))
        error('you must first install EasyH5 from http://github.com/fangq/easyh5/');
    end
    if(~isfield(opt,'rootname'))
        if(nargin>=3 && ischar(varargin{1}))
            opt.rootname=varargin{1};
        else
            opt.rootname=inputname(2);
        end
    end
    [varargout{1:nargout}]=saveh5(varargin{2},filename, opt);
else
    error('file suffix must be one of .json,.jnii,.jdt,.jmsh,.jnirs,.bjd,.bnii,.jdb,.bmsh,.bnirs,.ubj,.msgpack,.h5,.hdf5,.snirf');
end
