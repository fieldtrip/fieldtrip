% FT_PREAMBLE_INIT is a helper script that is used at the start of all FieldTrip main
% functions. It checks whether the user specified at lease one input arguments (i.e.
% the cfg) or shows the help of the calling function. It merges the global defaults
% with the cfg. It checks whether the output file already exists and whether it is OK
% to overwrite it. It tracks the function call.
%
% Use as
%   ft_preamble init
%
% See also FT_PREAMBLE, FT_POSTAMBLE, FT_TRACKUSAGE

% Copyright (C) 2011-2016, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% disabled for now, see further down
global ft_default

% Which high-level FieldTrip function is lowest on the stack, and which is
% highest? Lowest could e.g. be ft_selectdata called by another FieldTrip
% function, highest should always be the top-level FieldTrip function that
% was explicitly called by the user. These are used further down below,
% e.g. in the reproducescript functionality.
% (Note that we prepend random strings to local variables in pre- and
% postambles to avoid contaminating the caller's namespace.)
FjmoT6aA_stack = dbstack('-completenames');
% FjmoT6aA_stack(1) is this script
% FjmoT6aA_stack(2) is the calling ft_postamble function
% FjmoT6aA_stack(3) is the lowest FieldTrip function on the stack, the one currently being evaluated
% the highest FieldTrip function on the stack is given by the highest n
% such that FjmoT6aA_stack(n).path begins with the FieldTrip path
FjmoT6aA_lowest_ft = FjmoT6aA_stack(3).name;
FjmoT6aA_highest_ft = FjmoT6aA_stack(3).name;
[FjmoT6aA_ft_ver, FjmoT6aA_ft_path] = ft_version;
for FjmoT6aA_k = 3:numel(FjmoT6aA_stack)
  if startsWith(FjmoT6aA_stack(FjmoT6aA_k).file, FjmoT6aA_ft_path)
    FjmoT6aA_highest_ft = FjmoT6aA_stack(FjmoT6aA_k).name;
  else
    % we are operating under the assumption here that a private FieldTrip
    % function (i.e. one whose name does not start with "ft_") never calls
    % another high-level FieldTrip function.
    break;
  end
end
% check whether we are currently in a top-level user-called FT function
% (note that checking of names only probably is not sufficient, as recursion
% could occur)
FjmoT6aA_current_ft_toplevel = strcmp(FjmoT6aA_lowest_ft, FjmoT6aA_highest_ft) && ...
  (FjmoT6aA_k == 4 || numel(FjmoT6aA_stack) == 3);

if ft_nargin==0
  help(FjmoT6aA_lowest_ft);
  % throw the error as if it happened in the original function
  msg.message     = 'This function requires one or multiple input arguments, please refer to the documentation above';
  msg.identifier  = '';
  msg.stack       = FjmoT6aA_stack;
  error(msg);
end % if nargin

% check if there are fieldnames in the cfg that suggest as if the user
% erroneously inputted a data argument
checkdatafields = isfield(cfg, {'cfg' 'label' 'dimord' 'trialinfo' 'avg' 'powspctrm'});
if any(checkdatafields)
  help(FjmoT6aA_lowest_ft);
  % throw the error as if it happened in the original function
  msg.message     = 'It seems as if the first input argument is a FieldTrip data structure, while a cfg is expected';
  msg.identifier  = '';
  msg.stack       = FjmoT6aA_stack;
  error(msg);
end

% convert automatically from cell-array to structure
if iscell(cfg)
  cfg = ft_keyval2cfg(cfg);
end

% check that it is an struct or empty numeric array
assert(isstruct(cfg) || (isnumeric(cfg) && isempty(cfg)), 'The configuration must be a structure or empty');

% this script requires some options that can be user-specified, but otherwise are obtained from ft_default
% merge the default options into the configuration, except the preamble field which is used for passing arguments
cfg = mergeconfig(cfg, ft_default);

% determine whether function execution should be aborted or continued
if isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)
  assert(any(strcmp(fieldnames(cfg), 'outputfilepresent')), 'cfg.outputfilepresent is a required option, please see FT_DEFAULTS');
  % check whether the output file already exists
  if ischar(cfg.outputfile)
    chiL7fee_outputfile = {cfg.outputfile};
  else
    chiL7fee_outputfile = cfg.outputfile;
  end
  for i=1:numel(chiL7fee_outputfile)
    [p, f, x] = fileparts(chiL7fee_outputfile{i});
    if isempty(p)
      % the relative path was speciield
      chiL7fee_outputfile{i} = fullfile(pwd, chiL7fee_outputfile{i});
    else
      % the absolute path was specified
      chiL7fee_outputfile{i} = chiL7fee_outputfile{i};
    end
    if ~exist(chiL7fee_outputfile{i}, 'file')
      ft_abort = false;
    else
      % the output file exists, determine how to deal with it
      switch cfg.outputfilepresent
        case 'keep'
          if ft_nargout>0
            % continue executing the parent function
            ft_warning('output file %s is already present, but you also requested an output argument: continuing function execution', chiL7fee_outputfile{i});
            ft_abort = false;
          else
            % stop executing the parent function
            ft_warning('output file %s is already present: aborting function execution', chiL7fee_outputfile{i});
            ft_abort = true;
          end
        case 'overwrite'
          ft_warning('output file %s is already present: it will be overwritten', chiL7fee_outputfile{i});
          ft_abort = false;
        case 'error'
          ft_error('output file %s is already present', chiL7fee_outputfile{i});
        otherwise
          ft_error('invalid option for cfg.outputfilepresent');
      end % case
    end
  end % for each of the output files
else
  % there is no reason to abort execution
  ft_abort = false;
end % if chiL7fee_outputfile{i}

if isfield(cfg, 'reproducescript') && ~isempty(cfg.reproducescript)
  % the reproducescript code should only be executed in a top-level FT function
  if ~FjmoT6aA_current_ft_toplevel
    % we are in a FT function that was called by another FT function
    cfg = rmfield(cfg, 'reproducescript');
  else
    % we are in a top-level FT function
    if ~isfolder(cfg.reproducescript)
      mkdir(cfg.reproducescript);
    end
    % user-specified cfg.inputfile or cfg.outputfile are not compatible
    % with cfg.reproducescript functionality, throw error to make them
    % mutually exclusive
    if (isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)) || ...
       (isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile))
      ft_error('cfg.reproducescript cannot be used together with a user-specified cfg.inputfile or cfg.outputfile');
    end
    % this variable is used in loadvar, savevar and savefig
    Fief7bee_reproducescript = cfg.reproducescript;
    cfg = rmfield(cfg, 'reproducescript');
    % pause one second to ensure that subsequent file names (which contain the time stamp) are unique
    pause(1);
  end
end

if false
  % this is currently generating too much data and therefore disabled
  if isfield(cfg, 'trackusage') && ~(isequal(cfg.trackusage, false) || isequal(cfg.trackusage, 'no') || isequal(cfg.trackusage, 'off'))
    % track the usage of the calling function
    stack = dbstack('-completenames');
    % stack(1) is this script
    % stack(2) is the calling ft_postamble function
    % stack(3) is the main FieldTrip function that we are interested in
    ft_trackusage('function call', 'function', stack(3).name);
  end % if trackusage
end
