function res = time(this, ind, format)
% Method for getting the time axis
% FORMAT res = time(this, ind, format)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, Stefan Kiebel
% $Id: time.m 2949 2009-03-25 11:57:16Z vladimir $

if this.Nsamples>0
    res = (0:(this.Nsamples-1))./this.Fsample + this.timeOnset;
else
    res = [];
end

if exist('ind') == 1 && ~isempty(ind)
    res = res(ind);
end

if exist('format') ==1
    if strcmp(format, 'ms')
        res = res*1000;
    end
end