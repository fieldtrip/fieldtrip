function res = transformtype(this, newtype)
% Method for getting/setting type of transform
% FORMAT res = transformtype(this, name)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: transformtype.m 3350 2009-09-03 13:19:20Z vladimir $

if nargin == 1
    res = this.transform.ID;
else
    if strncmpi(this, 'TF', 2) && length(size(this))~=4
        error('TF transformtype can only be assigned to 4D dataset');
    end
    
    this.transform.ID = newtype;
    res = this;
end
