function res = repl(this, varargin)
% Method for getting replication counts, over trials
% FORMAT res = repl(this, index, nrepl)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: repl.m 2729 2009-02-11 10:37:25Z vladimir $

res = getset(this, 'trials', 'repl', varargin{:});