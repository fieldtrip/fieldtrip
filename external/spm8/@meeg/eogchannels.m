function ind = eogchannels(this)
% Return indices of EOG channels
% FORMAT ind = eogchannels(this)
%
%  this      - MEEG object
%  ind       - row vector of indices of EOG channels
%
% See also ecgchannels, emgchannels, meegchannels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: eogchannels.m 2884 2009-03-16 18:27:25Z guillaume $

type = chantype(this);
ind = find(ismember(upper(type), {'EOG', 'HEOG', 'VEOG'}));
ind = ind(:)'; % must be row to allow to use it as loop indices
