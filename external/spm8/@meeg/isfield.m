function res = isfield(this, varargin)
% Returns true if the string fieldname is the name of a field in the 
% substructure 'other' in the meeg object 'this'.
% FORMAT res = isfield(this,fieldname)
%
% An overloaded function...
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips
% $Id: isfield.m 2720 2009-02-09 19:50:46Z vladimir $

res = isfield(this.other, varargin{:});