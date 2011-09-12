function h = ft_connectivityplot(cfg, varargin)

% FT_CONNECTIVITYPLOT plots frequency domain connectivity with dimord
% 'chan_chan_freq'. The data are rendered in a square matrix of 'subplots'
% using fast ft_plot_vector for rendering.
%
% Use as
%   ft_connectivityplot(cfg, data)
%
% The data is a structure containing the output to FT_CONNECTIVITYANALYSIS
% using a frequency domain metric of connectivity, and with a dimord of
% 'chan_chan_freq'.
%
% The cfg can have the following options:
%   cfg.parameter   = string, the functional parameter to be plotted (default = 'cohspctrm')
%   cfg.zlim        = [lower upper]
%   cfg.channel     = list of channels to be included for the plotting (default = 'all')
%
% See also FT_CONNECTIVITYANALYSIS, FT_PLOT_VECTOR, FT_CHANNELSELECTION

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
% $Id: ft_connectivityplot$

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();;
ftFuncMem   = memtic();

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'renamed',	 {'zparam', 'parameter'});

cfg.channel   = ft_getopt(cfg, 'channel', 'all');
cfg.parameter = ft_getopt(cfg, 'parameter', 'cohspctrm');
cfg.zlim      = ft_getopt(cfg, 'zlim', []);
cfg.color     = ft_getopt(cfg, 'color', 'brgkywrgbkywrgbkywrgbkyw');

% make the function recursive if numel(varargin)>1
% FIXME check explicitly which channels belong together
if numel(varargin)>1
  data = varargin{1};
  tmpcfg = cfg;
  if ischar(cfg.parameter)
    % do nothing
  elseif iscell(cfg.parameter)
    tmpcfg.parameter = cfg.parameter{1};
  end
  ft_connectivityplot(tmpcfg, data);
  tmpcfg = cfg;
  for k = 2:numel(varargin)
    tmpcfg.color   = tmpcfg.color(2:end);
    tmpcfg.holdfig = 1;
    if ischar(cfg.parameter)
      % do nothing
    elseif iscell(cfg.parameter)
      tmpcfg.parameter = cfg.parameter{k};
    end
    ft_connectivityplot(tmpcfg, varargin{k});
  end
  return;
else
  data = varargin{1};
end

if strcmp(data.dimord, 'chan_chan_freq')
  % that's ok
elseif strcmp(data.dimord, 'chancmb_freq')
  % convert into 'chan_chan_freq'
  data = ft_checkdata(data, 'cmbrepresentation', 'full');
else
  error('the data should have a dimord of %s or %s', 'chan_chan_freq', 'chancmb_freq');
end

if ~isfield(data, cfg.parameter)
  error('the data does not contain the requested parameter %s', cfg.parameter);
end

cfg.channel = ft_channelselection(cfg.channel, data.label);
data        = ft_selectdata(data, 'channel', cfg.channel);

dat   = data.(cfg.parameter);
nchan = numel(data.label);
nfreq = numel(data.freq);

if isempty(cfg.zlim)
  cfg.zlim = [min(dat(:)) max(dat(:))];
end

if (isfield(cfg, 'holdfig') && cfg.holdfig==0) || ~isfield(cfg, 'holdfig')
  h = figure;hold on;
end

for k = 1:nchan
  for m = 1:nchan
    if k~=m
      ix  = k;
      iy  = nchan - m + 1;
      % use the convention of the row-channel causing the column-channel
      tmp = reshape(dat(m,k,:), [nfreq 1]);
      ft_plot_vector(tmp, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'vlim', cfg.zlim, 'box', 'yes', 'color', cfg.color(1));
    end
  end
end

for k = 1:nchan
  ft_plot_text(0,       (nchan + 1 - k).*1.2, data.label{k});
  ft_plot_text(k.*1.2,  (nchan + 1)    .*1.2, data.label{k});
end

axis([-0.2 (nchan+1).*1.2+0.2 0 (nchan+1).*1.2+0.2]);
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$Id: ft_movieplotER.m 4096 2011-09-03 15:49:40Z roboos $'; % this will be auto-updated by the revision control system

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();

% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername(); % this is helpful for debugging
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

