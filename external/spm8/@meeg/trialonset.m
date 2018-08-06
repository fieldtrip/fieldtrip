function res = trialonset(this, varargin)
% Method for getting/setting trial onset times
% FORMAT res = trialonset(this, ind, onset)
%   ind = indices of trials
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: trialonset.m 1373 2008-04-11 14:24:03Z spm $

res = getset(this, 'trials', 'onset', varargin{:});
