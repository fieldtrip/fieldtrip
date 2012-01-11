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

% determine the name of the calling FieldTrip function
s = dbstack;
fname = s(3).name;
clear s

% determine the input and output arguments to the calling FieldTrip function
[inarg, outarg] = funargname(fname);
fname   = str2func(fname);
argin   = cell(size(inarg));  % this will be filled further down
argout  = cell(size(outarg)); % this will be filled further down

% gather the input arguments
for i=1:length(inarg)
  eval(sprintf('argin{%d} = %s;', i, inarg{i}));
end
clear i

if strcmp(inarg(end), 'varargin')
  % the variable list of input arguments should not be passed as a single cell-array
  argin = cat(2, argin(1), argin{2});
end

for i=1:length(argin)
  % if one of the input arguments is being evaluated in the background,
  % then keep the subsequent processing also in the background
  if isa(argin{i}, 'background')
    cfg.distribute = 'background';
  end
end

if isfield(cfg, 'distribute') && ~isempty(cfg.distribute)
  
  % prevent re-distribution
  distribute = cfg.distribute;
  cfg = rmfield(cfg, 'distribute');
  
  % gather the input arguments for a second time, needed because the cfg has changed
  for i=1:length(inarg)
    eval(sprintf('argin{%d} = %s;', i, inarg{i}));
  end
  clear i

  if strcmp(inarg(end), 'varargin')
    % the variable list of input arguments should not be passed as a single cell-array
    argin = cat(2, argin(1), argin{2});
  end

  % get the options from the cfg as key-value pairs
  % these can be specified as sub-structure in cfg.peer, cfg.distcomp, etc.
  options = ft_getopt(cfg, distribute);
  options = ft_cfg2keyval(options);

  switch distribute
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'peer'
      % submit the job for remote execution
      jobid = peerfeval(fname, argin{:}, options{:});
      % collect the output arguments
      [argout{:}] = peerget(jobid, 'timeout', inf);
      
      % clean the local variables
      clear jobid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'qsub'
      % submit the job for remote execution
      jobid = qsubfeval(fname, argin{:}, options{:});
      % collect the output arguments
      [argout{:}] = qsubget(jobid, 'timeout', inf);
      
      % clean the local variables
      clear jobid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'distcomp'
      % each input arg should be a cell-array
      for i=1:length(argin)
        argin{i} = {argin{i}};
      end
      
      % submit the job for execution using the MATLAB distributed computing toolbox
      % and optionally MATLAB distributed computing engines
      jobid = dfevalasync(fname, length(argout), argin{:}, options{:});
      % collect the output arguments
      waitForState(jobid);
      [argout{:}] = getAllOutputArguments(jobid);
      
      % clean the local variables
      clear jobid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'background'
      % submit the job for background evaluation
      [argout{:}] = background(fname, argin{:}, options{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
      error('unsupported method "%s" for distributed computing', distribute);
      
  end % switch
  
  % reassign the output arguments
  for i=1:length(outarg)
    eval(sprintf('%s = argout{%d};', outarg{i}, i));
  end
  clear i
  
  % it is now safe to add the option back to the configuration
  cfg.distribute = distribute;
  
  % clean up the local variables
  clear fname distribute options
  
end % if distribute

clear inarg outarg argin argout
