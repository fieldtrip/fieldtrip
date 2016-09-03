function res = size(this, varargin)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 3767 2010-03-09 22:49:30Z vladimir $

res = size(this.data.y, varargin{:});

if ntrials(this) == 1 && isempty(varargin)
    res = [res 1];
end
