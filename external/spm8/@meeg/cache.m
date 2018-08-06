function res = cache(this, stuff)
% Method for retrieving/putting stuff from/into the temporary cache
% FORMAT res = cache(obj, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: cache.m 1373 2008-04-11 14:24:03Z spm $

if ischar(stuff)
    % assume this is a retrieval
        res = getcache(this, stuff);        
else
    % assume this is a put
    name = inputname(2);
        res = setcache(this, name, stuff);
end

function res = getcache(this, stuff)
if isfield(this.cache, stuff)
    eval(['res = this.cache.' stuff ';']);
else
    res = [];
end

function this = setcache(this, name, stuff)
eval(['this.cache(1).' name ' = stuff;']);
