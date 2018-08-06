function res = fnamedat(this, name)
% Method for getting/setting file name of data file
% FORMAT res = fnamedat(this, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: fnamedat.m 1373 2008-04-11 14:24:03Z spm $


switch nargin
    case 1
        res = getfnamedat(this);        
    case 2
        res = setfnamedat(this, name);
    otherwise
end

function res = getfnamedat(this)
res = this.data.fnamedat;

function this = setfnamedat(this, name)
this.data.fnamedat = name;
