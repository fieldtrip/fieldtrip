function debugCleanup

% DEBUGCLEANUP is a cleanup function that is being used by FT_PREAMBLE_DEBUG. It is
% called when a high-level FieldTrip function exits, either after finishing successfully or after detecting an error.
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
  
  case {'saveonerror' 'save' 'saveonsuccess'}
    fprintf('-----------------------------------------------------------------------\n');
    if strcmp(Ce9dei2ZOo_debug, 'saveonerror')
      fprintf('An error was detected while executing %s\n', Ce9dei2ZOo_funname);
    elseif strcmp(Ce9dei2ZOo_debug, 'saveonsuccess')
      fprintf('Execution of %s finished successfully\n', Ce9dei2ZOo_funname);
    end
    filename = fullfile(tempdir, [funname '_' datestr(now, 30) '.mat']);
    variables = cat(2, {'funname', 'last_err', 'last_error', 'last_warning'}, {Ce9dei2ZOo_argin.name});
    try
      fprintf('Saving debug information to %s\n', filename);
      save(filename, variables{:})
    catch
      % it might fail because the disk is full
      fprintf('error while saving debug information to %s\n', filename);
    end
    fprintf('-----------------------------------------------------------------------\n');
    
  case {'displayonerror' 'display' 'displayonsuccess'}
    fprintf('-----------------------------------------------------------------------\n');
    if strcmp(Ce9dei2ZOo_debug, 'displayonerror')
      fprintf('An error was detected while executing %s with\n', Ce9dei2ZOo_funname);
    elseif strcmp(Ce9dei2ZOo_debug, 'displayonsuccess')
      fprintf('Execution of %s finished successfully\n', Ce9dei2ZOo_funname);
    end
    fprintf('\n');
    for i=1:length(Ce9dei2ZOo_argin)
      fprintf('%s =\n', Ce9dei2ZOo_argin(i).name);
      switch class(Ce9dei2ZOo_argin(i).value)
        case 'config'
          disp(struct(Ce9dei2ZOo_argin(i).value));
        case 'char'
          fprintf('%s\n\n', Ce9dei2ZOo_argin(i).value); % with an extra newline
        otherwise
          disp(Ce9dei2ZOo_argin(i).value);
      end % switch class
    end
    fprintf('-----------------------------------------------------------------------\n');
    
  otherwise
    % do nothing
    
end % switch

% clean up the global variables
Ce9dei2ZOo_debug   = 'no';
Ce9dei2ZOo_funname = [];
Ce9dei2ZOo_argin   = [];
