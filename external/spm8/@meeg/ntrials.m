function res = ntrials(obj)
% Method for getting the number of trials in the file
% FORMAT res = ntrials(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: ntrials.m 1373 2008-04-11 14:24:03Z spm $

res = length(obj.trials);