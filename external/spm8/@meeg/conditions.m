function res = conditions(this, varargin)
% Method for getting condition labels, over trials
% FORMAT res = conditions(this, ind, conditionlabels)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: conditions.m 4310 2011-04-18 16:07:35Z guillaume $

res = getset(this, 'trials', 'label', varargin{:});

if nargin == 1 && ~iscell(res)
    res = {res};
end
