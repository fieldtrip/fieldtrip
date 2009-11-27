function [data] = om_load_full(filename,format)

% LOAD_BIN   Load full Matrix
%
%   Load full Matrix
%
%   SYNTAX
%       [DATA] = OM_LOAD_FULL(FILENAME,FORMAT)
%
%       FORMAT : can be 'ascii' or 'binary' (default)
%
%   Created by Alexandre Gramfort on 2007-11-27.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information

if nargin == 1
    format = 'binary';
end

switch format
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

