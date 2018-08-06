function [data] = om_load_sym(filename,format)

% OM_LOAD_SYM   Load symmetric Matrix
%
%   Load symmetric Matrix
%
%   SYNTAX
%       [DATA] = OM_LOAD_SYM(FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' or 'matlab' (default)
%

% Copyright (C) 2010-2017, OpenMEEG developers

me = 'OM_LOAD_SYM';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

if nargin == 1
    format = 'mat';
end

isOctave = exist('OCTAVE_VERSION') ~= 0;

switch format
    case 'matlab'
        if isOctave
            data_raw = load(filename);
        else
            data_raw = load(filename,'-mat');
        end
        if isfield(data_raw, 'symmatrix')
            sym = data_raw.symmatrix;
        elseif isfield(data_raw, 'data') && isfield(data_raw, 'size')
            sym = data_raw;
        end
        clear data_raw;
        if length(sym.data) ~= sym.size * (sym.size+1) / 2
            error('Number of entries in symmatrix doesn''t fit to the size of the matrix.');
        end
        data = zeros(sym.size);
        data(triu(ones(sym.size,sym.size)) > 0) = sym.data;
        data = data';
        data(triu(ones(sym.size,sym.size)) > 0) = sym.data;
        clear sym;
    case 'binary'
        file = fopen(filename,'r');
        dim = fread(file,1,'uint32','ieee-le');
        data = zeros(dim,dim);
        data(triu(ones(dim,dim)) > 0) = fread(file,dim*(dim+1)/2,'double','ieee-le');
        data = data + data' - diag(diag(data));
        fclose(file);
    case 'ascii'
        file = fopen(filename);
        rawdata = fscanf(file,'%f');
        % rawdata = cell2mat(rawdata);
        dim = (-1 + sqrt(1+8*length(rawdata)))/2;
        assert(dim == ceil(dim),'Bad dimension for a symmetric Matrix')
        data = zeros(dim,dim);
        data(tril(ones(dim,dim)) > 0) = rawdata;
        data = data + data' - diag(diag(data));
        fclose(file);
    otherwise
        error([me,' : Unknown file format'])
end

