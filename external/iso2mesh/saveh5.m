function saveh5(data, fname, varargin)
%
%    saveh5(data, outputfile)
%       or
%    saveh5(data, outputfile, options)
%    saveh5(data, outputfile, 'Param1',value1, 'Param2',value2,...)
%
%    Save a MATLAB struct (array) or cell (array) into an HDF5 file
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        data: a structure (array) or cell (array) to be stored.
%        fname: the output HDF5 (.h5) file name
%        options: (optional) a struct or Param/value pairs for user specified options
%            JData [0|1] use JData Specifiation to serialize complex data structures
%                         such as complex/sparse arrays, tables, maps, graphs etc by
%                         calling jdataencode before saving data to HDF5
%            RootName: the HDF5 path of the root object. If not given, the
%                         actual variable name for the data input will be used as
%                         the root object. The value shall not include '/'.
%            UnpackHex [1|0]: convert the 0x[hex code] in variable names
%                         back to Unicode string using decodevarname.m
%            Compression: ['deflate'|''] - use zlib-deflate method
%                         to compress data array
%            CompressArraySize: [100|int]: only to compress an array if the
%                         total element count is larger than this number.
%            CompressLevel: [5|int] - a number between 1-9 to set
%                         compression level
%            Chunk: a size vector or empty - breaking a large array into
%                         small chunks of size specified by this parameter
%            Append [0|1]: if set to 1, new data will be appended to an
%                         file if already exists under the user-defined
%                         'rootname' path; if set to 0, old data
%                         will be overwritten if the file exists.
%            VariableLengthString [0|1]: if set to 1, strings and char arrays will be
%                         saved with type 'H5T_C_S1' and size 'H5T_VARIABLE'
%            Scalar [1|0]: if set to 1, arrays of length 1 will be saved as
%                         a scalar instead of a length-1 array
%            Transpose: [1|0] - if set to 1 (default), MATLAB arrays are
%                         transposed (from column-major to row-major) so
%                         that the output HDF5 dataset shows the same
%                         dimensions as in MATLAB when reading from other
%                         tools.
%            ComplexFormat: {'realKey','imagKey'}: use 'realKey' and 'imagKey'
%                  as keywords for the real and the imaginary part of a
%                  complex array, respectively (sparse arrays not supported);
%                  the default values are {'Real','Imag'}
%
%    example:
%        a=struct('a',rand(5),'b','string','c',true,'d',2+3i,'e',{'test',[],1:5});
%        saveh5(a,'test.h5');
%        saveh5(a(1),'test2.h5','rootname','');
%        saveh5(a(1),'test2.h5','compression','deflate','compressarraysize',1);
%        saveh5(a,'test.h5j','jdata',1);
%        saveh5(a,'test.h5j','rootname','/newroot','append',1);
%
%    this file is part of EasyH5 Toolbox: https://github.com/NeuroJSON/easyh5
%
%    License: GPLv3 or 3-clause BSD license, see https://github.com/NeuroJSON/easyh5 for details
%

if (nargin < 2)
    error('you must provide at least two inputs');
end

rootname = ['/' inputname(1)];
opt = struct;

if (length(varargin) == 1 && ischar(varargin{1}))
    rootname = [varargin{1} '/' inputname(1)];
else
    opt = varargin2struct(varargin{:});
end

opt.compression = jsonopt('Compression', '', opt);
opt.compresslevel = jsonopt('CompressLevel', 5, opt);
opt.compressarraysize = jsonopt('CompressArraySize', 100, opt);
opt.unpackhex = jsonopt('UnpackHex', 1, opt);
opt.dotranspose = jsonopt('Transpose', 1, opt);
opt.variablelengthstring = jsonopt('VariableLengthString', 0, opt);
opt.scalar = jsonopt('Scalar', 1, opt);

opt.releaseid = 0;
vers = ver('MATLAB');
if (~isempty(vers))
    opt.releaseid = datenum(vers(1).Date);
end
opt.skipempty = (opt.releaseid < datenum('1-Jan-2015'));

if (isfield(opt, 'rootname'))
    rootname = ['/' opt.rootname];
end

if (~isfield(opt, 'rootname') && ~isempty(regexp(rootname, '/$', 'once')))
    rootname = [rootname 'data'];
end

if (jsonopt('JData', 0, opt))
    data = jdataencode(data, 'Base64', 0, 'UseArrayZipSize', 0, opt);
end

if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
    save(fname, 'data', '-hdf5');
    return
end

try
    if (isa(fname, 'H5ML.id'))
        fid = fname;
    else
        if (jsonopt('append', 0, opt) && exist(fname, 'file'))
            fid = H5F.open(fname, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
        else
            fid = H5F.create(fname, 'H5F_ACC_TRUNC', H5P.create('H5P_FILE_CREATE'), H5P.create('H5P_FILE_ACCESS'));
        end
    end
    obj2h5(rootname, data, fid, 1, opt);
catch ME
    if (exist('fid', 'var') && fid > 0)
        H5F.close(fid);
    end
    rethrow(ME);
end

if (~isa(fname, 'H5ML.id'))
    H5F.close(fid);
end

%% -------------------------------------------------------------------------
function oid = obj2h5(name, item, handle, level, varargin)

if (iscell(item))
    oid = cell2h5(name, item, handle, level, varargin{:});
elseif (isstruct(item))
    oid = struct2h5(name, item, handle, level, varargin{:});
elseif (ischar(item) || isa(item, 'string'))
    oid = mat2h5(name, item, handle, level, varargin{:});
elseif (isa(item, 'containers.Map'))
    oid = map2h5(name, item, handle, level, varargin{:});
elseif (isa(item, 'categorical'))
    oid = cell2h5(name, cellstr(item), handle, level, varargin{:});
elseif (islogical(item) || isnumeric(item) || isa(item, 'timeseries'))
    oid = mat2h5(name, item, handle, level, varargin{:});
else
    oid = any2h5(name, item, handle, level, varargin{:});
end

%% -------------------------------------------------------------------------
function oid = idxobj2h5(name, idx, varargin)
oid = obj2h5(sprintf('%s%d', name, idx), varargin{:});

%% -------------------------------------------------------------------------
function oid = cell2h5(name, item, handle, level, varargin)

num = numel(item);
if (num > 1)
    idx = reshape(1:num, size(item));
    idx = num2cell(idx);
    oid = cellfun(@(x, id) idxobj2h5(name, id, x, handle, level, varargin{:}), item, idx, 'UniformOutput', false);
else
    oid = cellfun(@(x) obj2h5(name, x, handle, level, varargin{:}), item, 'UniformOutput', false);
end

%% -------------------------------------------------------------------------
function oid = struct2h5(name, item, handle, level, varargin)

num = numel(item);
if (num > 1)
    oid = obj2h5(name, num2cell(item), handle, level, varargin{:});
else
    pd = 'H5P_DEFAULT';
    gcpl = H5P.create('H5P_GROUP_CREATE');
    tracked = H5ML.get_constant_value('H5P_CRT_ORDER_TRACKED');
    indexed = H5ML.get_constant_value('H5P_CRT_ORDER_INDEXED');
    order = bitor(tracked, indexed);
    H5P.set_link_creation_order(gcpl, order);
    if (varargin{1}.unpackhex)
        name = decodevarname(name);
    end
    try
        handle = H5G.create(handle, name, pd, gcpl, pd);
        isnew = 1;
    catch
        isnew = 0;
    end

    names = fieldnames(item);
    oid = cell(1, length(names));
    for i = 1:length(names)
        oid{i} = obj2h5(names{i}, item.(names{i}), handle, level + 1, varargin{:});
    end

    if (isnew)
        H5G.close(handle);
    end
end

%% -------------------------------------------------------------------------
function oid = map2h5(name, item, handle, level, varargin)

pd = 'H5P_DEFAULT';
gcpl = H5P.create('H5P_GROUP_CREATE');
tracked = H5ML.get_constant_value('H5P_CRT_ORDER_TRACKED');
indexed = H5ML.get_constant_value('H5P_CRT_ORDER_INDEXED');
order = bitor(tracked, indexed);
H5P.set_link_creation_order(gcpl, order);
try
    if (varargin{1}.unpackhex)
        name = decodevarname(name);
    end
    handle = H5G.create(handle, name, pd, gcpl, pd);
    isnew = 1;
catch
    isnew = 0;
end

names = item.keys;
oid = zeros(length(names));
for i = 1:length(names)
    oid(i) = obj2h5(names{i}, item(names{i}), handle, level + 1, varargin{:});
end

if (isnew)
    H5G.close(handle);
end

%% -------------------------------------------------------------------------
function oid = mat2h5(name, item, handle, level, varargin)
typemap = h5types;

opt = varargin{1};

is1dvector = 0;

if (isa(item, 'timeseries'))
    if (item.TimeInfo.Length == 1 || (item.TimeInfo.isUniform && item.TimeInfo.Increment == 1 && ndims(item.Data) == 3 && size(item.Data, 1) == 1 && size(item.Data, 2) == 1))
        is1dvector = 1;
        item = squeeze(item.Data);
    else
        item = [item.Time, item.Data];
    end
end

if (opt.dotranspose)
    item = permute(item, ndims(item):-1:1);
end

pd = 'H5P_DEFAULT';
gcpl = H5P.create('H5P_GROUP_CREATE');
tracked = H5ML.get_constant_value('H5P_CRT_ORDER_TRACKED');
indexed = H5ML.get_constant_value('H5P_CRT_ORDER_INDEXED');
order = bitor(tracked, indexed);
H5P.set_link_creation_order(gcpl, order);

if (~(isfield(opt, 'complexformat') && iscellstr(opt.complexformat) && numel(opt.complexformat) == 2) || strcmp(opt.complexformat{1}, opt.complexformat{2}))
    opt.complexformat = {'Real', 'Imag'};
end

usefilter = opt.compression;
complevel = opt.compresslevel;
minsize = opt.compressarraysize;
chunksize = jsonopt('Chunk', size(item), opt);

if (isa(item, 'logical'))
    item = uint8(item);
end

if (~isempty(usefilter) && numel(item) >= minsize)
    if (isnumeric(usefilter) && usefilter(1) == 1)
        usefilter = 'deflate';
    end
    if (strcmpi(usefilter, 'deflate'))
        pd = H5P.create('H5P_DATASET_CREATE');
        h5_chunk_dims = fliplr(chunksize);
        H5P.set_chunk(pd, h5_chunk_dims);
        H5P.set_deflate(pd, complevel);
        opt.scalar = 0;
        opt.variablelengthstring = 0;
    else
        error('Filter %s is unsupported', usefilter);
    end
end

if (opt.unpackhex)
    name = decodevarname(name);
end

oid = [];

if (isempty(item) && opt.skipempty)
    warning('The HDF5 library is older than v1.8.7, and can not save empty datasets. Skip saving "%s"', name);
    return
end

if (isreal(item) || isa(item, 'string'))
    if (issparse(item))
        idx = find(item);
        oid = sparse2h5(name, struct('Size', size(item), 'SparseIndex', idx, 'Real', full(item(idx))), handle, level, varargin{:});
    else
        itemtype = H5T.copy(typemap.(class(item)));
        if ((ischar(item) || isa(item, 'string')) && opt.variablelengthstring)
            H5T.set_size(itemtype, 'H5T_VARIABLE');
            itemsize = H5S.create('H5S_SCALAR');
            if (isa(item, 'string') && is1dvector)
                itemsize = H5S.create_simple(1, length(item), length(item));
            elseif (isa(item, 'string') && length(item) > 1)
                itemsize = H5S.create_simple(ndims(item), fliplr(size(item)), fliplr(size(item)));
            end
        elseif (isnumeric(item) && numel(item) == 1 && ndims(item) == 2 && opt.scalar && ~is1dvector)
            itemsize = H5S.create('H5S_SCALAR');
        else
            if (is1dvector)
                itemsize = H5S.create_simple(1, length(item), length(item));
            else
                itemsize = H5S.create_simple(ndims(item), fliplr(size(item)), fliplr(size(item)));
            end
        end
        try
            oid = H5D.create(handle, name, itemtype, itemsize, pd);
        catch ME
            if (jsonopt('append', 0, opt) == 1)
                H5L.delete(handle, name, 'H5P_DEFAULT');
                oid = H5D.create(handle, name, itemtype, H5S.create_simple(ndims(item), fliplr(size(item)), fliplr(size(item))), pd);
            end
        end
        if ((ischar(item) || isa(item, 'string')) && isfield(opt, 'variablelengthstring')  && opt.variablelengthstring)
            if (isa(item, 'string'))
                H5D.write(oid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', cellstr(item));
            else
                H5D.write(oid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', {item});
            end
        else
            H5D.write(oid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', item);
        end
    end
else
    if (issparse(item))
        idx = find(item);
        oid = sparse2h5(name, struct('Size', size(item), 'SparseIndex', idx, 'Real', full(real(item(idx))), 'Imag', full(imag(item(idx)))), handle, level, varargin{:});
    else
        typeid = H5T.copy(typemap.(class(item)));
        elemsize = H5T.get_size(typeid);
        memtype = H5T.create ('H5T_COMPOUND', elemsize * 2);
        H5T.insert (memtype, opt.complexformat{1}, 0, typeid);
        H5T.insert (memtype, opt.complexformat{2}, elemsize, typeid);
        if (is1dvector)
            oid = H5D.create(handle, name, memtype, H5S.create_simple(1, length(item), length(item)), pd);
        else
            oid = H5D.create(handle, name, memtype, H5S.create_simple(ndims(item), fliplr(size(item)), fliplr(size(item))), pd);
        end
        H5D.write(oid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', struct(opt.complexformat{1}, real(item), opt.complexformat{2}, imag(item)));
    end
end
if (~isempty(oid))
    H5D.close(oid);
end

%% -------------------------------------------------------------------------
function oid = sparse2h5(name, item, handle, level, varargin)

opt = varargin{1};

idx = item.SparseIndex;

if (isempty(idx) && opt.skipempty)
    warning('The HDF5 library is older than v1.8.7, and can not save empty datasets. Skip saving "%s"', name);
    oid = [];
    return
end

adata = item.Size;
item = rmfield(item, 'Size');
hasimag = isfield(item, 'Imag');

typemap = h5types;

pd = 'H5P_DEFAULT';

usefilter = opt.compression;
complevel = opt.compresslevel;
minsize = opt.compressarraysize;
chunksize = jsonopt('Chunk', size(item), opt);

if (~isempty(usefilter) && numel(idx) >= minsize)
    if (isnumeric(usefilter) && usefilter(1) == 1)
        usefilter = 'deflate';
    end
    if (strcmpi(usefilter, 'deflate'))
        pd = H5P.create('H5P_DATASET_CREATE');
        h5_chunk_dims = fliplr(chunksize);
        H5P.set_chunk(pd, h5_chunk_dims);
        H5P.set_deflate(pd, complevel);
    else
        error('Filter %s is unsupported', usefilter);
    end
end

idxtypeid = H5T.copy(typemap.(class(idx)));
idxelemsize = H5T.get_size(idxtypeid);
datatypeid = H5T.copy(typemap.(class(item.Real)));
dataelemsize = H5T.get_size(datatypeid);
memtype = H5T.create ('H5T_COMPOUND', idxelemsize + dataelemsize * (1 + hasimag));
H5T.insert (memtype, 'SparseIndex', 0, idxtypeid);
H5T.insert (memtype, 'Real', idxelemsize, datatypeid);
if (hasimag)
    H5T.insert (memtype, 'Imag', idxelemsize + dataelemsize, datatypeid);
end
oid = H5D.create(handle, name, memtype, H5S.create_simple(ndims(idx), fliplr(size(idx)), fliplr(size(idx))), pd);
H5D.write(oid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', item);

space_id = H5S.create_simple(ndims(adata), fliplr(size(adata)), fliplr(size(adata)));
attr_size = H5A.create(oid, 'SparseArraySize', H5T.copy('H5T_NATIVE_DOUBLE'), space_id, H5P.create('H5P_ATTRIBUTE_CREATE'));
H5A.write(attr_size, 'H5ML_DEFAULT', adata);
H5A.close(attr_size);

%% -------------------------------------------------------------------------
function oid = any2h5(name, item, handle, level, varargin)
pd = 'H5P_DEFAULT';

if (varargin{1}.unpackhex)
    name = decodevarname(name);
end

rawdata = getByteStreamFromArray(item);  % use undocumented matlab function
oid = H5D.create(handle, name, H5T.copy('H5T_STD_U8LE'), H5S.create_simple(ndims(rawdata), size(rawdata), size(rawdata)), pd);
H5D.write(oid, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', pd, rawdata);

adata = class(item);
space_id = H5S.create_simple(ndims(adata), size(adata), size(adata));
attr_type = H5A.create(oid, 'MATLABObjectClass', H5T.copy('H5T_C_S1'), space_id, H5P.create('H5P_ATTRIBUTE_CREATE'));
H5A.write(attr_type, 'H5ML_DEFAULT', adata);
H5A.close(attr_type);

adata = size(item);
space_id = H5S.create_simple(ndims(adata), size(adata), size(adata));
attr_size = H5A.create(oid, 'MATLABObjectSize', H5T.copy('H5T_NATIVE_DOUBLE'), space_id, H5P.create('H5P_ATTRIBUTE_CREATE'));
H5A.write(attr_size, 'H5ML_DEFAULT', adata);
H5A.close(attr_size);

H5D.close(oid);

%% -------------------------------------------------------------------------
function typemap = h5types
typemap.char = 'H5T_C_S1';
typemap.string = 'H5T_C_S1';
typemap.double = 'H5T_IEEE_F64LE';
typemap.single = 'H5T_IEEE_F32LE';
typemap.logical = 'H5T_STD_U8LE';
typemap.uint8 = 'H5T_STD_U8LE';
typemap.int8 = 'H5T_STD_I8LE';
typemap.uint16 = 'H5T_STD_U16LE';
typemap.int16 = 'H5T_STD_I16LE';
typemap.uint32 = 'H5T_STD_U32LE';
typemap.int32 = 'H5T_STD_I32LE';
typemap.uint64 = 'H5T_STD_U64LE';
typemap.int64 = 'H5T_STD_I64LE';
