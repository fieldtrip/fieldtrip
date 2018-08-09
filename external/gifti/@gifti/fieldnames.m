function names = fieldnames(this)
% Fieldnames method for GIfTI objects
% FORMAT names = fieldnames(this)
% this   -  GIfTI object
% names  -  field names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: fieldnames.m 6507 2015-07-24 16:48:02Z guillaume $

if numel(this) > 1, warning('Only handle scalar objects yet.'); end

pfn = {'vertices','faces','normals','cdata','mat','labels','indices'};

names = unique(pfn(isintent(this,pfn)));
