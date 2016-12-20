function res = fsample(this, value)
% Method for getting and setting the sampling rate
% FORMAT res = fsample(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fsample.m 3200 2009-06-12 17:29:40Z vladimir $

if nargin == 1
    res = this.Fsample;
else
    this.Fsample = value;
    res = this;
end