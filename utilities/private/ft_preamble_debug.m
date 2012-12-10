% FT_PREAMBLE_DEBUG is a helper script for debugging problems with FieldTrip functions
%
% In case of an error, the debugCleanup function does the remainder of the work.
%
% In case a normal termination, pt_postamble_debug will be called first, followed by
% debugCleanup.
%
% See also FT_POSTAMBLE_DEBUG, DEBUGCLEANUP

% these variables are shared by the three debug handlers
global Ce9dei2ZOo_debug Ce9dei2ZOo_funname Ce9dei2ZOo_argin 

Ce9dei2ZOo_debug   = 'no';
Ce9dei2ZOo_funname = [];
Ce9dei2ZOo_argin   = [];

if ~isfield(cfg, 'debug')
  % do not provide extra debugging facilities
  return
end

% reset the last error and warning
lasterr('');
lastwarn('');

% remember the input arguments
Ce9dei2ZOo_workspace = evalin('caller', 'whos');
Ce9dei2ZOo_workspace = setdiff({Ce9dei2ZOo_workspace.name}, {'ft_default', 'revision', 'ftohDiW7th_FuncClock', 'ftohDiW7th_FuncMem', 'ftohDiW7th_FuncTimer', 'Ce9dei2ZOo_debug', 'Ce9dei2ZOo_funname', 'Ce9dei2ZOo_argin'});
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
    Ce9dei2ZOo_debug = 'none';
    
  case 'saveonerror'
    Ce9dei2ZOo_debug = 'saveonerror';
    onCleanup(@debugCleanup);
    
  case 'display'
    Ce9dei2ZOo_debug = 'display';
    debugCleanup; % call it once
    Ce9dei2ZOo_debug = 'none';
    
  case 'displayonerror'
    Ce9dei2ZOo_debug = 'displayonerror';
    onCleanup(@debugCleanup);
    
  otherwise
    % do nothing
    
end % switch
