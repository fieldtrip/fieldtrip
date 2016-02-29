function [stat] = ft_networkanalysis(cfg, data)

% FT_NETWORKANALYSIS computes various network graph measures from
% between-channel or between source-level EEG/MEG signals. This function
% acts as a wrapper aroun the network metrics implemented in the brain
% connectivity toolbox developed by Olaf Sporns and colleagues.
%
% Use as
%   stat = ft_networkanalysis(cfg, data)
%
% where the first input argument is a configuration structure (see below)
% and the second argument is the output of FT_CONNECTIVITYANALYSIS.
%
% At present the input data should be channel-level data with dimord
% 'chan_chan(_freq)(_time)' or source data with dimord
% 'pos_pos(_freq)(_time)'.
%
% The configuration structure has to contain
%   cfg.method    = string, specifying the graph measure that will be
%                   computed. See below for the list of supported measures.
%   cfg.parameter = string specifying the bivariate parameter in the data
%                   for which the graph measure will be computed.
%
% Supported methods are
%   assortativity
%   betweenness,      betweenness centrality (nodes)
%   charpath,         characteristic path length, needs distance matrix as
%                     input
%   clustering_coef,  clustering coefficient
%   degrees
%   density
%   distance
%   edge_betweenness, betweenness centrality (edges)
%   transitivity
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% 	cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a
% *.mat file on disk and/or the output data will be written to a *.mat
% file. These mat files should contain only a single variable,
% corresponding with the input/output structure.
%
% See also FT_CONNECTIVITYANALYSIS, FT_CONNECTIVITYPLOT

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

cfg = ft_checkconfig(cfg, 'required', {'method' 'parameter'});

cfg.threshold = ft_getopt(cfg, 'threshold', []);

% ensure that the bct-toolbox is on the path
ft_hastoolbox('BCT', 1);

% check the data for the correct dimord and for the presence of the requested parameter
if ~strcmp(data.dimord(1:7), 'pos_pos') && ~strcmp(data.dimord(1:9), 'chan_chan'),
  error('the dimord of the input data should start with ''chan_chan'' or ''pos_pos''');
end

% conversion to double is needed because some BCT functions want to do matrix
% multiplications on boolean matrices
input = double(data.(cfg.parameter));

% some metrics explicitly require a certain parameter
if strcmp(cfg.method, 'charpath') && ~strcmp(cfg.parameter, 'distance')
  error('characteristic path length can only be computed on distance matrices');
end

% check for binary or not
isbinary = true;
for k = 1:size(input,3)
  for m = 1:size(input,4)
    tmp = input(:,:,k,m);
    isbinary = all(ismember(tmp(:), [0 1]));
    if ~isbinary,
      break;
    end
  end
end
if isbinary
  fprintf('input graph is binary\n');
else
  fprintf('input graph is weighted\n');
end

if ~isbinary && ~isempty(cfg.threshold)
  fprintf('thresholding the input graph at a value of %d\n', cfg.threshold);
  newinput = false(size(input));
  for k = 1:size(input,3)
    for m = 1:size(input,4)
      tmp = input(:,:,k,m);
      newinput(:,:,k,m) = tmp>cfg.threshold;
    end
  end
  input = double(newinput); clear newinput;
  isbinary = true;
end

% check for directed or not
isdirected = true;
for k = 1:size(input,3)
  for m = 1:size(input,4)
    tmp = input(:,:,k,m);
    isdirected = ~all(all(tmp==tmp.'));
    if ~isdirected,
      break;
    end
  end
end
if isdirected
  fprintf('input graph is directed\n');
else
  fprintf('input graph is undirected\n');
end

fprintf('computing %s\n', cfg.method);
% allocate memory
needlabel = true;
switch cfg.method
  case {'assortativity' 'charpath' 'density'  'transitivity'}
    % 1 value per connection matrix
    outsiz = [size(input) 1];
    outsiz(1:2) = [];
    output = zeros(outsiz);
    if strcmp(data.dimord(1:3), 'pos')
      dimord = data.dimord(9:end);
    elseif strcmp(data.dimord(1:4), 'chan')
      dimord = data.dimord(11:end);
    end
    needlabel = false;
  case {'betweenness' 'clustering_coef' 'degrees'}
    % 1 value per node
    outsiz = [size(input) 1];
    outsiz(1) = [];
    output = zeros(outsiz);
    if strcmp(data.dimord(1:3), 'pos')
      dimord = data.dimord(5:end);
    elseif strcmp(data.dimord(1:4), 'chan')
      dimord = data.dimord(6:end);
    end
  case {'distance' 'edge_betweenness'}
    % 1 value per node pair
    outsiz = [size(input) 1];
    output = zeros(outsiz);
    dimord = data.dimord;
end

binarywarning = 'weights are not taken into account and graph is converted to binary values by thresholding';

for k = 1:size(input, 3)
  for m = 1:size(input, 4)

    % switch to the appropriate function from the BCT
    switch cfg.method
      case 'assortativity'
        if ~isbinary, ft_warning(binarywarning); end

        if isdirected
          output(k,m) = assortativity(input(:,:,k,m), 1);
        elseif ~isdirected
          output(k,m) = assortativity(input(:,:,k,m), 0);
        end
      case 'betweenness'
        if isbinary
          output(:,k,m) = betweenness_bin(input(:,:,k,m));
        elseif ~isbinary
          output(:,k,m) = betweenness_wei(input(:,:,k,m));
        end
      case 'breadthdist'
        error('not yet implemented');
      case 'charpath'
        % this needs the distance matrix as input, this is dealt with
        % above
        output(:,k) = charpath(input(:,:,k,m))';
      case 'clustering_coef'
        if isbinary && isdirected
          output(:,k,m) = clustering_coef_bd(input(:,:,k,m));
        elseif isbinary && ~isdirected
          output(:,k,m) = clustering_coef_bu(input(:,:,k,m));
        elseif ~isbinary && isdirected
          output(:,k,m) = clustering_coef_wd(input(:,:,k,m));
        elseif ~isbinary && ~isdirected
          output(:,k,m) = clustering_coef_wu(input(:,:,k,m));
        end
      case 'degrees'
        if ~isbinary, ft_warning(binarywarning); end

        if isdirected
          [in, out, output(:,k,m)] = degrees_dir(input(:,:,k,m));
          % FIXME do something here
        elseif ~isdirected
          output(:,k,m) = degrees_und(input(:,:,k,m));
        end
      case 'density'
        if ~isbinary, ft_warning(binarywarning); end

        if isdirected
          output(k,m) = density_dir(input(:,:,k,m));
        elseif ~isdirected
          output(k,m) = density_und(input(:,:,k,m));
        end
      case 'distance'
        if isbinary
          output(:,:,k,m) = distance_bin(input(:,:,k,m));
        elseif ~isbinary
          output(:,:,k,m) = distance_wei(input(:,:,k,m));
        end
      case 'edge_betweenness'
        if isbinary
          output(:,:,k,m) = edge_betweenness_bin(input(:,:,k,m));
        elseif ~isbinary
          output(:,:,k,m) = edge_betweenness_wei(input(:,:,k,m));
        end
      case 'efficiency'
        error('not yet implemented');
      case 'modularity'
        error('not yet implemented');
      case 'participation_coef'
        error('not yet implemented');
      case 'transitivity'
        if isbinary && isdirected
          output(k,m) = transitivity_bd(input(:,:,k,m));
        elseif isbinary && ~isdirected
          output(k,m) = transitivity_bu(input(:,:,k,m));
        elseif ~isbinary && isdirected
          output(k,m) = transitivity_wd(input(:,:,k,m));
        elseif ~isbinary && ~isdirected
          output(k,m) = transitivity_wu(input(:,:,k,m));
        end
      otherwise
        error('unsupported connectivity metric %s requested');
    end

  end % for m
end % for k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the output structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stat              = [];
stat.(cfg.method) = output;
stat.dimord       = dimord;
if isfield(data, 'label') && needlabel,  stat.label  = data.label;  end
if isfield(data, 'freq'),   stat.freq   = data.freq;   end
if isfield(data, 'time'),   stat.time   = data.time;   end
if isfield(data, 'grad'),   stat.grad   = data.grad;   end
if isfield(data, 'elec'),   stat.elec   = data.elec;   end
if exist('dof',  'var'),    stat.dof    = dof;         end
% FIXME this needs to be implemented still

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance stat
ft_postamble history    stat
ft_postamble savevar    stat
