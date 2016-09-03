function res = path(this, name)
% Method for getting/setting path
% FORMAT res = path(this, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: path.m 1373 2008-04-11 14:24:03Z spm $

switch nargin
    case 1
        res = getpath(this);        
    case 2
        res = setpath(this, name);
    otherwise
end

function res = getpath(this)
try
    res = this.path;
catch
    res = '';
end

function this = setpath(this, name)
this.path = name;


