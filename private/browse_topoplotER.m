function browse_topoplotER(cfg, data)

% this is a simple helper function for DATABROWSER

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: browse_topoplotER.m,v $
% Revision 1.1  2009/10/19 14:19:21  roboos
% first version, to work with databrowser
%

% convert to an ERP
timelock = timelockanalysis([], data);

default             = [];
default.xlim        = [min(timelock.time) max(timelock.time)];
default.marker      = 'on';
default.interactive = 'no';
cfg = mergeconfig(cfg, default);

figure; 
topoplotER(cfg, timelock);
