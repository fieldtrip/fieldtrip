function res = getcache(obj, name)
% Method for retrieving stuff from the temporary cache
% FORMAT res = getcache(obj, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: getcache.m 1373 2008-04-11 14:24:03Z spm $

if isfield(obj.cache, name)
    eval(['res = obj.cache.' name ';']);
else
    res = [];
end
