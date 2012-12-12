function debugCleanup

% DEBUGCLEANUP is a cleanup function that is being used by FT_PREAMBLE_DEBUG. It is
% called when a high-level FieldTrip function exits, either normally or in error.
%
% See also FT_PREAMBLE_DEBUG, FT_POSTAMBLE_DEBUG

% these variables are shared by the three debug handlers
global Ce9dei2ZOo_debug Ce9dei2ZOo_funname Ce9dei2ZOo_argin

funname = Ce9dei2ZOo_funname;
for i=1:length(Ce9dei2ZOo_argin)
  eval(sprintf('%s = Ce9dei2ZOo_argin(%d).value;', Ce9dei2ZOo_argin(i).name, i));
end
last_err     = lasterr;
last_error   = lasterror;
last_warning = lastwarn;

switch Ce9dei2ZOo_debug
  
  case {'saveonerror' 'save'}
    filename = fullfile(tempdir, [funname '_' datestr(now, 30) '.mat']);
    variables = cat(2, {'funname', 'last_err', 'last_error', 'last_warning'}, {Ce9dei2ZOo_argin.name});
    try
      save(filename, variables{:})
      fprintf('saving debug information to %s\n', filename);
    catch
      % it might fail because the disk is full
      fprintf('error while saving debug information to %s\n', filename);
    end
    
  case {'displayonerror' 'display'}
    fprintf('\n');
    for i=1:length(Ce9dei2ZOo_argin)
      fprintf('%s =\n', Ce9dei2ZOo_argin(i).name);
      disp(Ce9dei2ZOo_argin(i).value);
    end
    
  otherwise
    % do nothing
    
end % switch

% clean up the global variables
Ce9dei2ZOo_debug   = 'no';
Ce9dei2ZOo_funname = [];
Ce9dei2ZOo_argin   = [];
