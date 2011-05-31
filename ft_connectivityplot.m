function h = ft_connectivityplot(cfg, data)

% ft_connectivityplot plots frequency domain connectivity with dimord
% 'chan_chan_freq'. The data are rendered in a square matrix of 'subplots'
% using fast ft_plot_vector for rendering.
%
% Use as:
%   ft_connectivityplot(cfg, data)
%
% The data is a structure containing the output to FT_CONNECTIVITYANALYSIS
% using a frequency domain metric of connectivity, and with a dimord of
% 'chan_chan_freq'.
%
% The cfg can have the following options:
%   cfg.zparam  = string (default = 'cohspctrm'), the functional parameter
%                  to be plotted.
%   cfg.channel = list of channels to be included for the plotting (default
%                   = 'all');
%
% See also:
%   FT_CONNECTIVITYANALYSIS, FT_PLOT_VECTOR, FT_CHANNELSELECTION

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

cfg.channel = ft_getopt(cfg, 'channel', 'all');
cfg.zparam  = ft_getopt(cfg, 'zparam', 'cohspctrm');

if ~strcmp(data.dimord, 'chan_chan_freq')
  error('the data should have a dimord of %s', 'chan_chan_freq');
end

if ~isfield(data, cfg.zparam)
  error('the data does not contain the requested zparam %s', cfg.zparam);
end

cfg.channel = ft_channelselection(cfg.channel, data.label);
data        = ft_selectdata(data, 'channel', cfg.channel);

dat   = data.(cfg.zparam);
nchan = numel(data.label);
nfreq = numel(data.freq);

h = figure;hold on;
for k = 1:nchan
  for m = 1:nchan
    if k~=m
      ix  = k;
      iy  = nchan - m + 1;
      tmp = reshape(dat(k,m,:), [nfreq 1]);
      ft_plot_vector(tmp, 'width', 1, 'height', 1, 'hpos', ix.*1.2, 'vpos', iy.*1.2, 'vlim', cfg.zlim, 'box', 'yes');
    end
  end
end

for k = 1:nchan
  ft_plot_text(0, (nchan - k + 1).*1.2, data.label{k});
  ft_plot_text(k.*1.2,  (nchan + 1).*1.2, data.label{k});
end

axis equal;
axis off;