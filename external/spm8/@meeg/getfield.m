function res = getfield(this, varargin)
% Returns  fields in .other
% FORMAT res = getfield(this, varargin)
%
% An overloaded function...
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: getfield.m 3228 2009-06-26 17:43:19Z vladimir $

res = getfield(this.other, varargin{:});