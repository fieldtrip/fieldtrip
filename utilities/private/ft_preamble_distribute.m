% FT_PREAMBLE_DISTRIBUTE is a helper script for distributed computing on a
% cluster using the MATLAB distributed computing toolbox, the FieldTrip peer
% toolbox or the FieldTrip qsub toolbox.

% Copyright (C) 2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

return; %JM added this because of unexpected behavior in FT-functions with optional input arguments (that are not defined

% determine the name of the calling FieldTrip function
s = dbstack;
fname = s(3).name;
clear s

% determine the input and output arguments to the calling FieldTrip function
[arginname, argoutname] = funargname(fname);
fname   = str2func(fname);
arginval   = cell(size(arginname));  % this will be filled further down
argoutval  = cell(size(argoutname)); % this will be filled further down

% it might be that some input arguments are optional
arginname = arginname(1:nargin);
arginval  = arginval(1:nargin);

% gather the input arguments
for i=1:length(arginname)
  eval(sprintf('arginval{%d} = %s;', i, arginname{i}));
end
clear i

if strcmp(arginname(end), 'varargin')
  % the variable list of input arguments should not be passed as a single cell-array
  arginval = cat(2, arginval(1), arginval{2});
end

for i=1:length(arginval)
  % if one of the input arguments is being evaluated in the background,
  % then keep the subsequent processing also in the background
  if isa(arginval{i}, 'background')
    cfg.distribute = 'background';
  end
end

if isfield(cfg, 'distribute') && ~isempty(cfg.distribute)
  
  % prevent re-distribution
  distribute = cfg.distribute;
  cfg = rmfield(cfg, 'distribute');
  
  % gather the input arguments for a second time, needed because the cfg has changed
  for i=1:length(arginname)
    eval(sprintf('arginval{%d} = %s;', i, arginname{i}));
  end
  clear i

  if strcmp(arginname(end), 'vararginval')
    % the variable list of input arguments should not be passed as a single cell-array
    arginval = cat(2, arginval(1), arginval{2});
  end

  % get the options from the cfg as key-value pairs
  % these can be specified as sub-structure in cfg.peer, cfg.distcomp, etc.
  options = ft_getopt(cfg, distribute);
  options = ft_cfg2keyval(options);

  switch distribute
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'local'
      % submit the job for background evaluation
      [argoutval{:}] = feval(fname, arginval{:});
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'peer'
      % submit the job for remote execution
      jobid = peerfeval(fname, arginval{:}, options{:});
      % collect the output arguments
      [argoutval{:}] = peerget(jobid, 'timeout', inf);
      
      % clean the local variables
      clear jobid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'qsub'
      % submit the job for remote execution
      jobid = qsubfeval(fname, arginval{:}, options{:});
      % collect the output arguments
      [argoutval{:}] = qsubget(jobid, 'timeout', inf);
      
      % clean the local variables
      clear jobid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'distcomp'
      % each input arg should be a cell-array
      for i=1:length(arginval)
        arginval{i} = {arginval{i}};
      end
      
      % submit the job for execution using the MATLAB distributed computing toolbox
      % and optionally MATLAB distributed computing engines
      jobid = dfevalasync(fname, length(argoutval), arginval{:}, options{:});
      % collect the output arguments
      waitForState(jobid);
      [argoutval{:}] = getAllOutputArguments(jobid);
      
      % clean the local variables
      clear jobid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'background'
      % submit the job for background evaluation
      [argoutval{:}] = background(fname, arginval{:}, options{:});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
      error('unsupported method "%s" for distributed computing', distribute);
      
  end % switch
  
  % reassign the output arguments
  for i=1:length(argoutname)
    eval(sprintf('%s = argoutval{%d};', argoutname{i}, i));
  end
  clear i
  
  % it is now safe to add the option back to the configuration
  cfg.distribute = distribute;
  
  % clean up the local variables
  clear fname distribute options
  
end % if distribute

clear arginname argoutname arginval argoutval
