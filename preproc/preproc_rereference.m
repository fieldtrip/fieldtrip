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
% $Log: preproc_rereference.m,v $
% Revision 1.3  2009/01/07 12:44:18  roboos
% also allow refchan='all'
%
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%
% Revision 1.1  2008/05/23 06:54:22  roboos
% created initial scratch version of preprocessing module, to be used in fieldtrip or as stand-alone toolbox (e.g. in spm8 or braingain)
% some functions are copies of existing roboos/misc versions, some just contain some example code for the implementation
%
% Revision 1.2  2003/03/17 10:37:28  roberto
% improved general help comments and added copyrights
%

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
