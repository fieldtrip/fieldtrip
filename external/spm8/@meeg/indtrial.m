function res = indtrial(this, label, flag)
% Method for getting trial indices based on condition labels
% FORMAT res = indtrial(this, label)
% this       - MEEG object
% label      - string or cell array of labels, 'GOOD' and 'BAD'
%              can be added to list of labels to select only
%              good or bad trials respectively
% flag       - 'GOOD' or 'BAD' to include only good or bad trials
%              respectively (all are selected by default)
%
% res        - vector of trial indices matching condition labels
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: indtrial.m 5059 2012-11-15 13:48:35Z vladimir $

if ischar(label)
    label = {label};
end

[dum, res] = match_str(label, conditions(this));

if nargin > 2
    if strcmpi(flag, 'GOOD')
        res = setdiff(res, badtrials(this));
    elseif strcmpi(flag, 'BAD')
        res = intersect(res, badtrials(this));
    end
end