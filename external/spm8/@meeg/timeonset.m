function res = timeonset(this, value)
% Method for reading and setting the time onset
% FORMAT res = timeonset(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: timeonset.m 2116 2008-09-18 14:50:34Z stefan $

if nargin == 1
    res = this.timeOnset;
else
    this.timeOnset = value;
    res = this;
end