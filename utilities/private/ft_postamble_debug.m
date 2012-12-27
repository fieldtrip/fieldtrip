% ft_postamble_debug is a helper script for debugging problems with fieldtrip functions
%
% see also ft_preamble_debug debugcleanup

% these variables are shared by the three debug handlers
global ce9dei2zoo_debug ce9dei2zoo_funname ce9dei2zoo_argin

if ~isfield(cfg, 'debug')
  % do not provide extra debugging facilities
  return
end

switch cfg.debug
  case {'displayonsuccess' 'saveonsuccess'}
    % do not clean up the global variables yet
    % these are still needed by the cleanup function
    
  otherwise
    % clean up the global variables
    % this results in the cleanup function doing nothing
    ce9dei2zoo_debug   = 'no';
    ce9dei2zoo_funname = [];
    ce9dei2zoo_argin   = [];
end
