% FT_PREAMBLE_DEBUG is a helper script for debugging problems with FieldTrip
% functions. It installs an "onCleanup" function that gets executed if an error is
% detected, allowing to automatically save the input data and configuration to disk.
%
% In case of an error, the DEBUGCLEANUP function does the remainder of the work. In
% case a normal termination, FT_POSTAMBLE_DEBUG will be called first, followed by
% DEBUGCLEANUP.
%
% Use as
%   ft_preamble debug
%   .... regular code goes here ...
%   ft_postamble debug
%
% See also FT_POSTAMBLE_DEBUG, DEBUGCLEANUP

% these variables are shared by the three debug handlers
global Ce9dei2ZOo_debug Ce9dei2ZOo_funname Ce9dei2ZOo_argin Ce9dei2ZOo_ws Ce9dei2ZOo_ns Ce9dei2ZOo_is Ce9dei2ZOo_ds

if ~isempty(Ce9dei2ZOo_debug) && ~isequal(Ce9dei2ZOo_debug, 'no')
  % the debugging handler is already set by a higher-level function
  return
end

if isfield(cfg, 'verbose') && ischar(cfg.verbose)
  % store the current state of the notifications
  Ce9dei2ZOo_ws = ft_warning;
  Ce9dei2ZOo_ns = ft_notice;
  Ce9dei2ZOo_is = ft_info;
  Ce9dei2ZOo_ds = ft_debug;
  switch cfg.verbose
    case 'error'
      ft_warning off
      ft_notice  off
      ft_info    off
      ft_debug   off
    case 'warning'
      ft_warning on
      ft_notice  off
      ft_info    off
      ft_debug   off
    case 'notice'
      ft_warning on
      ft_notice  on
      ft_info    off
      ft_debug   off
    case 'info'
      ft_warning on
      ft_notice  on
      ft_info    on
      ft_debug   off
    case 'debug'
      ft_warning on
      ft_notice  on
      ft_info    on
      ft_debug   on
    otherwise
      % just leave them as they are
  end
end

if ~isfield(cfg, 'debug')
  % do not provide extra debugging facilities
  return
end

% reset the global variables used to handle the debugging
Ce9dei2ZOo_debug   = 'no';
Ce9dei2ZOo_funname = [];
Ce9dei2ZOo_argin   = [];

% reset the last error and warning
lasterr('');
lastwarn('');

% remember the variables that were passed as input arguments
Ce9dei2ZOo_workspace = evalin('caller', 'whos');
Ce9dei2ZOo_workspace = Ce9dei2ZOo_workspace(~strcmp({Ce9dei2ZOo_workspace.class}, 'function_handle')); % only variables, not anonymous functions
Ce9dei2ZOo_workspace = setdiff({Ce9dei2ZOo_workspace.name}, {'ft_default', 'ft_revision', 'ft_nargin', 'ft_abort', 'ftohDiW7th_FuncMem', 'ftohDiW7th_FuncTimer', 'Ce9dei2ZOo_debug', 'Ce9dei2ZOo_funname', 'Ce9dei2ZOo_argin'});
Ce9dei2ZOo_argin     = [];
for Ce9dei2ZOo_indx=1:length(Ce9dei2ZOo_workspace)
  Ce9dei2ZOo_argin(Ce9dei2ZOo_indx).name  = Ce9dei2ZOo_workspace{Ce9dei2ZOo_indx};
  Ce9dei2ZOo_argin(Ce9dei2ZOo_indx).value = evalin('caller', Ce9dei2ZOo_workspace{Ce9dei2ZOo_indx});
end
clear Ce9dei2ZOo_workspace Ce9dei2ZOo_indx

% stack(1) is this script
% stack(2) is the calling ft_postamble function
% stack(3) is the main FieldTrip function that we are interested in
Ce9dei2ZOo_funname = dbstack;
Ce9dei2ZOo_funname = Ce9dei2ZOo_funname(3).name;

switch cfg.debug
  case 'save'
    Ce9dei2ZOo_debug = 'save';
    debugCleanup; % call it once
    Ce9dei2ZOo_debug = 'no';
    
  case 'saveonerror'
    Ce9dei2ZOo_debug  = 'saveonerror';
    % when the Ce9dei2ZOo_handle gets deleted, the debugCleanup function will be executed
    Ce9dei2ZOo_handle = onCleanup(@debugCleanup);
    
  case 'saveonsuccess'
    Ce9dei2ZOo_debug  = 'saveonsuccess';
    % when the Ce9dei2ZOo_handle gets deleted, the debugCleanup function will be executed
    Ce9dei2ZOo_handle = onCleanup(@debugCleanup);
    
  case 'display'
    Ce9dei2ZOo_debug = 'display';
    debugCleanup; % call it once
    Ce9dei2ZOo_debug = 'no';
    
  case 'displayonerror'
    Ce9dei2ZOo_debug  = 'displayonerror';
    % when the Ce9dei2ZOo_handle gets deleted, the debugCleanup function will be executed
    Ce9dei2ZOo_handle = onCleanup(@debugCleanup);
    
  case 'displayonsuccess'
    Ce9dei2ZOo_debug  = 'displayonsuccess';
    % when the Ce9dei2ZOo_handle gets deleted, the debugCleanup function will be executed
    Ce9dei2ZOo_handle = onCleanup(@debugCleanup);
    
  otherwise
    % do nothing
    
end % switch
