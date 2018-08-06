function res = pickconditions(this, label, rejectbad)
% Method for returning indices of trials of a certain trial type.
% If input argument rejectbad == 0, the function will also return trial 
% indices which are set to bad (i.e. rejected). If bad is omitted, 
% the default is not to include rejected trials.
% FORMAT res = pickconditions(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: pickconditions.m 3146 2009-05-26 09:54:23Z vladimir $

if nargin<3
    rejectbad = 1;
end

c = conditions(this);

if isa(label, 'char')
    label = {label};
end

res = find(ismember(c, label));

if rejectbad && ~isempty(res)
    res = res(~reject(this, res));
end