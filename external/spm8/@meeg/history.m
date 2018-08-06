function res = history(this, varargin)
% Method for getting or adding to the history of function calls of some
% M/EEG data
% FORMAT res = history(this, varargin)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: history.m 3196 2009-06-11 12:54:47Z vladimir $


if isempty(varargin)
    res = this.history;
else
    % add another history item
    if length(varargin) > 2 % To enable reset
        nh = 0;
        this.history = [];
    else
        nh = length(this.history);
    end

    if ischar(varargin{1})
        this.history(nh+1).fun = varargin{1};

        if isstruct(varargin{2}) && isfield(varargin{2}, 'D') && ...
                isa(varargin{2}.D, 'meeg')
            varargin{2}.D = fullfile(varargin{2}.D.path, varargin{2}.D.fname);
        end

        this.history(nh+1).args = varargin{2};

    elseif isstruct(varargin{1})
        this.history = varargin{1};
    end

    res = this;
end