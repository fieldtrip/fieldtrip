function res = units(this, varargin)
% Method for setting/getting all units, over channels
% FORMAT res = units(this, ind)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: units.m 1373 2008-04-11 14:24:03Z spm $

res = getset(this, 'channels', 'units', varargin{:});
