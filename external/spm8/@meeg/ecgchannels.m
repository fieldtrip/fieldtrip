function ind = ecgchannels(this)
% Return indices of ECG channels
% FORMAT ind = ecgchannels(this)
%
%  this      - MEEG object
%  ind       - row vector of indices of ECG channels
%
% See also eogchannels, emgchannels, meegchannels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: ecgchannels.m 2884 2009-03-16 18:27:25Z guillaume $

type = chantype(this);
ind = find(ismember(upper(type), {'ECG', 'EKG'}));
ind = ind(:)'; % must be row to allow to use it as loop indices
