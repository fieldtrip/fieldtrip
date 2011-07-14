function watchdog

% WATCHDOG will trigger an exit() if the master disappears or if the allowed time elapsed
%
% To enable the watchdog you should do 
%   watchdog(masterid, time)
% and to disable it again
%   clear watchdog

warning('could not locate mex file');

