function res = fieldnames(this, varargin)
% Returns names of the fields in .other
% FORMAT res = fieldnames(this)
%
% An overloaded function...
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fieldnames.m 2720 2009-02-09 19:50:46Z vladimir $

res = fieldnames(this.other, varargin{:});