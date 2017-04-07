function res = reject(this, varargin)
% Method for getting/setting rejection flags
% FORMAT res = reject(this, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: reject.m 1681 2008-05-19 12:32:08Z vladimir $


res = getset(this, 'trials', 'bad', varargin{:});

if iscell(res)
    res = [res{:}];
end

