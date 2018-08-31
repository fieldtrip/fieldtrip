% FT_PREAMBLE_INIT is a helper script that is used at the start of all FieldTrip main
% functions. It checks whether the user specified at lease one input arguments (i.e.
% the cfg) or shows the help of the calling function. It merges the global defaults
% with the cfg. It checks whether the output file already exists and whether it is OK
% to overwrite it. It tracks the function call.
%
% Use as
%   ft_preamble init
%
% See also FT_TRACKUSAGE

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

if ft_nargin==0
  stack = dbstack('-completenames');
  % stack(1) is this script
  % stack(2) is the calling ft_postamble function
  % stack(3) is the main FieldTrip function that we are interested in
  stack = stack(3);
  help(stack.name);
  % throw the error as if it happened in the original function
  msg.message     = 'This function requires one or multiple input arguments, please refer to the documentation above';
  msg.identifier  = '';
  msg.stack       = stack;
  error(msg);
end % if nargin

% check if there are fieldnames in the cfg that suggest as if the user
% erroneously inputted a data argument
checkdatafields = isfield(cfg, {'cfg' 'label' 'dimord' 'trialinfo' 'avg' 'powspctrm'});
if any(checkdatafields)
  stack = dbstack('-completenames');
  stack = stack(3);
  help(stack.name);
  % throw the error as if it happened in the original function
  msg.message     = 'It seems as if the first input argument is a FieldTrip data structure, while a cfg is expected';
  msg.identifier  = '';
  msg.stack       = stack;
  error(msg);
end

% convert automatically from cell-array to structure
if iscell(cfg)
  cfg = ft_keyval2cfg(cfg);
end

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
          ft_error('invalid option for cfg.chiL7fee_outputfile{i}present');
      end % case
    end
  end % for each of the output files
else
  % there is no reason to abort execution
  ft_abort = false;
end % if chiL7fee_outputfile{i}

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
