function tf = isfield(this,field)
% Isfield method for GIfTI objects
% FORMAT tf = isfield(this,field)
% this   -  GIfTI object
% field  -  string of cell array
% tf     -  logical array
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: isfield.m 6507 2015-07-24 16:48:02Z guillaume $

tf = ismember(field, fieldnames(this));
