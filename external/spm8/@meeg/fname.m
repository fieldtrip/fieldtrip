function res = fname(this, name)
% Method for getting/setting file name
% FORMAT res = fname(this, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: fname.m 1373 2008-04-11 14:24:03Z spm $

switch nargin
    case 1
        res = getfname(this);        
    case 2
        res = setfname(this, name);
    otherwise
end

function res = getfname(this)
res = this.fname;

function this = setfname(this, name)
this.fname = name;
