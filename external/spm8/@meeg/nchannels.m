function res = nchannels(this)
% returns number of channels
% FORMAT res = nchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nchannels.m 1373 2008-04-11 14:24:03Z spm $

res = length(this.channels);
