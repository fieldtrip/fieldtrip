function this = rmfield(this, fields)
% Method for removing an object field
% FORMAT this = rmfield(this, fields)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: rmfield.m 2883 2009-03-16 11:58:48Z vladimir $

 this.other = rmfield(this.other, fields);
