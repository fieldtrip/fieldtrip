function res = events(this, varargin)
% Method for getting/setting events per trial
% FORMAT res = events(this, ind, event)
%   ind = indices of trials
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: events.m 1373 2008-04-11 14:24:03Z spm $

res = getset(this, 'trials', 'events', varargin{:});
