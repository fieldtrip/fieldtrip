function [data] = om_load_sym(filename,format)

% OM_LOAD_SYM   Load symmetric Matrix
%
%   Load symmetric Matrix
%
%   SYNTAX
%       [DATA] = OM_LOAD_SYM(FILENAME,FORMAT)
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
    dim = fread(file,1,'uint32','ieee-le');
    data = zeros(dim,dim);
    data(triu(ones(dim,dim)) > 0) = fread(file,dim*(dim+1)/2,'double','ieee-le');
    data = data + data' - diag(diag(data));
    fclose(file);
case 'ascii'
    file = fopen(filename);
    rawdata = textscan(file,'%f');
    rawdata = cell2mat(rawdata);
    dim = (-1 + sqrt(1+8*length(rawdata)))/2;
    assert(dim == ceil(dim),'Bad dimension for a symmetric Matrix')
    data = zeros(dim,dim);
    data(tril(ones(dim,dim)) > 0) = rawdata;
    data = data + data' - diag(diag(data));
    fclose(file);
otherwise
    error([me,' : Unknown file format'])
end

