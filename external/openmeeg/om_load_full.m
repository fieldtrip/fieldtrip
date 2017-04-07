function [data] = om_load_full(filename,format)

% LOAD_BIN   Load full Matrix
%
%   Load full Matrix
%
%   SYNTAX
%       [DATA] = OM_LOAD_FULL(FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' or 'matlab' (default)
%

% Copyright (C) 2010-2017, OpenMEEG developers

me = 'OM_LOAD_FULL';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

if nargin == 1
    format = 'matlab';
end

isOctave = exist('OCTAVE_VERSION') ~= 0;

switch format
    case 'matlab'
        if isOctave
            data_raw = load(filename);
        else
            data_raw = load(filename,'-mat');
        end
        data = data_raw.linop;
        clear data_raw;
    case 'binary'
        file = fopen(filename,'r');
        dims = fread(file,2,'uint32','ieee-le');
        data = fread(file,prod(dims),'double','ieee-le');
        data = reshape(data,dims');
        fclose(file);
    case 'ascii'
        data = load(filename);
    otherwise
        error([me,' : Unknown file format'])
end
