function names = fieldnames(this)
% Fieldnames method for GIfTI objects
% FORMAT names = fieldnames(this)
% this   -  GIfTI object
% names  -  field names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$

if numel(this) > 1, warning('Only handle scalar objects yet.'); end

pfn = {'vertices' 'faces' 'normals' 'cdata' 'mat' 'labels'};

names = pfn(isintent(this,pfn));
