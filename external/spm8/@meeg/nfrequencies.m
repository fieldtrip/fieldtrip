function res = nfrequencies(this)
% Method for getting the number of frequencies for TF data
% FORMAT res = nsamples(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: nfrequencies.m 2846 2009-03-10 17:38:32Z guillaume $

if ~strncmp(transformtype(this), 'TF',2)
    res = [];
else
    res = length(this.transform.frequencies);
end