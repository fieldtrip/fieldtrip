function browse_multiplotER(cfg, data)

% this is a simple helper function for DATABROWSER

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: browse_multiplotER.m,v $
% Revision 1.1  2009/10/19 14:19:21  roboos
% first version, to work with databrowser
%

% convert to an ERP
timelock = timelockanalysis([], data);

default             = [];
default.interactive = 'yes';
cfg = mergeconfig(cfg, default);

figure;
multiplotER(cfg, timelock);
