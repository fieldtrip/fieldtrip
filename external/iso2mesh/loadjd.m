function varargout = loadjd(filename, varargin)
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
%                *.json,.jnii,.jdt,.jmsh,.jnirs,.jbids: JSON/JData based data files, see https://neurojson.org/jdata/draft2
%                *.bjd,.bnii,.jdb,.bmsh,.bnirs,.pmat: binary JData (BJData) files, see https://neurojson.org/bjdata/draft2
%                *.ubj: UBJSON-encoded files, see http://ubjson.org
%                *.msgpack: MessagePack-encoded files, see http://msgpack.org
%                *.h5,.hdf5,.snirf,.nwb: HDF5 files, see https://www.hdfgroup.org/
%                *.nii,.nii.gz: NIfTI files, need http://github.com/NeuroJSON/jnifty
%                *.tsv,.tsv.gz,.csv,.csv.gz: TSV/CSV files, need http://github.com/NeuroJSON/jbids
%                *.bval,.bvec: EEG .bval and .bvec files
%                *.mat: MATLAB/Octave .mat files
%        options: (optional) for JSON/JData files, these are optional 'param',value pairs
%                supported by loadjson.m; for BJData/UBJSON/MessagePack files, these are
%                options supported by loadbj.m; for HDF5 files, these are options
%                supported by loadh5.m (part of EasyH5 toolbox, http://github.com/NeuroJSON/easyh5/)
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

if (nargin < 1)
    error('you must provide file name');
end

if (regexpi(filename, '\.json$|\.jnii$|\.jdt$|\.jdat$|\.jmsh$|\.jnirs|\.jbids$'))
    [varargout{1:nargout}] = loadjson(filename, varargin{:});
elseif (regexpi(filename, '\.bjd$|\.bnii$|\.jdb$|\.jbat$|\.bmsh$|\.bnirs$|\.pmat$'))
    [varargout{1:nargout}] = loadbj(filename, varargin{:});
elseif (regexpi(filename, '\.ubj$'))
    [varargout{1:nargout}] = loadubjson(filename, varargin{:});
elseif (regexpi(filename, '\.msgpack$'))
    [varargout{1:nargout}] = loadmsgpack(filename, varargin{:});
elseif (regexpi(filename, '\.h5$|\.hdf5$|\.snirf$|\.nwb$'))
    if (~exist('loadh5', 'file'))
        error('you must first install EasyH5 from http://github.com/NeuroJSON/easyh5/');
    end
    [varargout{1:nargout}] = loadh5(filename, varargin{:});
elseif (regexpi(filename, '\.nii$|\.nii\.gz$'))
    if (~exist('loadnifti', 'file'))
        error('you must first install JNIFTY toolbox from http://github.com/NeuroJSON/jnifty/');
    end
    [varargout{1:nargout}] = loadnifti(filename, varargin{:});
elseif (regexpi(filename, '\.tsv$|\.tsv\.gz$|\.csv$|\.csv\.gz$'))
    if (~exist('loadbidstsv', 'file'))
        error('you must first install JBIDS toolbox from http://github.com/NeuroJSON/jbids/');
    end
    delim = sprintf('\t');
    if (regexpi(filename, '\.csv'))
        delim = ',';
    end
    [varargout{1:nargout}] = loadbidstsv(filename, delim);
elseif (regexpi(filename, '\.mat$|\.bvec$|\.bval$'))
    [varargout{1:nargout}] = load(filename, varargin{:});
else
    warning('only support parsing .json,.jnii,.jdt,.jmsh,.jnirs,.jbids,.bjd,.bnii,.jdb,.bmsh,.bnirs,.ubj,.msgpack,.h5,.hdf5,.snirf,.pmat,.nwb,.nii,.nii.gz,.tsv,.tsv.gz,.csv,.csv.gz,.mat,.bvec,.bval; load unparsed raw data');
    [varargout{1:nargout}] = fileread(filename);
end
