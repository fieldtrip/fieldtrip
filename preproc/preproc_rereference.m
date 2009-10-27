function [dat, ref] = preproc_rereference(dat, refchan)

% PREPROC_REREFERENCE computes the average reference over all EEG channels
% or rereferences the data to the selected channeld
%
% Use as
%  	[dat] = preproc_rereference(dat, refchan)
% where
%   dat        data matrix (Nchans X Ntime)
%   refchan    vector with indices of the new reference channels
%
% If the new reference channel is not specified, the data will be
% rereferenced to the average of all channels.
%
% See also PREPROC

% Copyright (C) 1998-2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% determine the new reference channels
if nargin<2 || isempty(refchan) || (ischar(refchan) && strcmp(refchan, 'all'))
  refchan = 1:Nchans;
end

% compute the average value over the reference channels
ref = mean(dat(refchan,:), 1);

% apply the new reference to the data
for chan=1:Nchans
  dat(chan,:) = dat(chan,:) - ref;
end
