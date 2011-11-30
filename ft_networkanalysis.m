function [stat] = ft_networkanalysis(cfg, data)

% FT_NETWORKANALYIS computes various graph measures from 
% between MEG/EEG channels or between source-level signals.
%
% Use as
%   stat = ft_networkanalysis(cfg, data)
% where the first input argument is a configuration structure (see
% below) and the second argument is the output of FT_CONNECTIVITYANALYSIS. 
% At present the input data should be channel level data with dimord 
% 'chan_chan(_freq)(_time)'. The metrics are computed using the brain 
% connectivity toolbox developed by Olaf Sporns and
% colleagues: www.brain-connectivity-toolbox.net.
%
% The configuration structure has to contain
%   cfg.method  = 
%
%     'clustering_coef'  clustering coefficient.
%     'degrees'
%
% Additional configuration options are
%   cfg.parameter   = string specifying the metric
% 
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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
% $Id: ft_connectivityanalysis.m 4812 2011-11-25 13:32:07Z jansch $

revision = '$Id: ft_connectivityanalysis.m 4812 2011-11-25 13:32:07Z jansch $';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

cfg = ft_checkconfig(cfg, 'required', 'method');

% set the defaults
cfg.inputfile   = ft_getopt(cfg, 'inputfile',   []);
cfg.outputfile  = ft_getopt(cfg, 'outputfile',  []);
cfg.parameter   = ft_getopt(cfg, 'parameter',   []);

ft_hastoolbox('BCT', 1);
if ~strcmp(data.dimord(1:9), 'chan_chan'),
  error('the dimord of the input data should start with ''chan_chan''');
end
if isempty(cfg.parameter)
  error('when using functions from the brain connectivity toolbox you should give a parameter for which the metric is to be computed');
else
  inparam = cfg.parameter;
end
outparam = cfg.method;

% gateway function to the brain connectivity toolbox
[datout, outdimord] = ft_connectivity_bct(data.(inparam), 'method', cfg.method, 'dimord', data.dimord);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the output structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stat            = [];
stat.(outparam) = datout;
stat.dimord     = outdimord;
if isfield(data, 'label'),  stat.label  = data.label;  end
if isfield(data, 'freq'),   stat.freq   = data.freq;   end
if isfield(data, 'time'),   stat.time   = data.time;   end
if isfield(data, 'grad'),   stat.grad   = data.grad;   end
if isfield(data, 'elec'),   stat.elec   = data.elec;   end
if exist('dof',  'var'),    stat.dof    = dof;         end
% FIXME this needs to be implemented still

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data
ft_postamble history stat
ft_postamble savevar stat
