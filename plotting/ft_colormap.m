function cmap = ft_colormap(varargin)

% FT_COLORMAP is a wrapper function with the same usage as the normal COLORMAP
% function, but it also knows about the colormaps from BREWERMAP and some colormaps
% from MATPLOTLIB. The recommended colormaps include 'parula', 'cividis', 'balance',
% and '*RdBu'.
%
% Use as
%   ft_colormap(name)
%   ft_colormap(name, n)
%   ft_colormap(handle, name)
%   ft_colormap(handle, name, n)
%
% The name is a string that specifies the colormap (see below). The optional handle
% can be used to specify the current figure (which is the default, see GCF) or the
% current axes (see GCA). The optional parameter n determines the number of steps or
% unique colors in the map (by default 64).
%
% The colormaps from MATLAB include 'parula', 'jet', 'hsv', 'hot', 'cool', 'spring',
% 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines',
% 'colorcube', 'prism', and 'flag'.
%
% The colormaps from MATPLOTLIB include 'cividis', 'inferno', 'magma', 'plasma',
% 'tab10', 'tab20', 'tab20b', 'tab20c', 'twilight', and 'viridis'.
%
% The colormaps from BREWERMAP include 'BrBG', 'PRGn', 'PiYG', 'PuOr', 'RdBu',
% 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1',
% 'Pastel2', 'Set1', 'Set2', 'Set3', 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens',
% 'Greys', 'OrRd', 'Oranges', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
% 'YlGn', 'YlGnBu', 'YlOrBr', and 'YlOrRd', plus their reverse when prefixed with '*'.
%
% The colormaps from CMOCEAN include 'thermal', 'haline', 'solar', 'ice', 'gray',
% 'oxy', 'deep', 'dense', 'algae', 'matter', 'turbid', 'speed', 'amp', 'tempo',
% 'rain', 'phase', 'topo', 'balance', 'delta', 'curl', 'diff', and 'tarn'.
%
% The colormaps from COLORCET include 'blueternary', 'coolwarm', 'cyclicgrey',
% 'depth', 'divbjy', 'fire', 'geographic', 'geographic2', 'gouldian', 'gray',
% 'greenternary', 'grey', 'heat', 'phase2', 'phase4', 'rainbow', 'rainbow2',
% 'rainbow3', 'rainbow4', 'redternary', 'reducedgrey', 'yellowheat', and all the ones
% with symbolic names.
%
% To reverse any of these these colormaps you can add a minus sign in front, like
% '-phase', '-balance' or '-RdBu'.
%
% Relevant publications:
% - Crameri et al. 2020. The misuse of colour in science communication. https://doi.org/10.1038/s41467-020-19160-7
% - Cooper et al. 2021. Over the rainbow: Guidelines for meaningful use of colour maps in neurophysiology. https://doi.org/10.1016/j.neuroimage.2021.118628
% - Kovesi 2015, Good colour maps: How to design them. https://doi.org/10.48550/arXiv.1509.03700
%
% See also COLORMAP, COLORMAPEDITOR, BREWERMAP, MATPLOTLIB, CMOCEAN, COLORCET

% Copyright (C) 2022, Jan-Mathijs Schoffelen, Robert Oostenveld
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

ft_hastoolbox('brewermap', 2);  % check and give a warning if it cannot be added
ft_hastoolbox('matplotlib', 2); % check and give a warning if it cannot be added
ft_hastoolbox('cmocean', 2);
ft_hastoolbox('colorcet', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpret the input arguments

% please note that ishandle('parula') returns a logical array
% hence the all(ishandle(...)), which in that case returns false

if nargin==0
  % apply the default colormap
  handle = gcf;
  name   = 'default';
  n      = 64;

elseif nargin==1 && isnumeric(varargin{1})
  % the user specified an Nx3 array as colormap
  handle = gcf;
  name   = varargin{1}; % note that this is not a string, it is dealt with below
  n      = nan;

elseif nargin==1 && all(ishandle(varargin{1}))
  % apply the default colormap on the specified handle
  handle = varargin{1};
  name   = 'default';
  n      = 64;

elseif nargin==1 && ischar(varargin{1})
  % apply the specified colormap on the current figure
  handle = gcf;
  name   = varargin{1};
  n      = 64;

elseif nargin==2 && all(ishandle(varargin{1}))
  % apply the specified colormap on the specified handle
  handle = varargin{1};
  name   = varargin{2};
  n      = 64;

elseif nargin==2 && ischar(varargin{1})
  % apply the specified colormap with specified N on the current figure
  handle = gcf;
  name   = varargin{1};
  n      = varargin{2};

elseif nargin==3 && all(ishandle(varargin{1}))
  % apply the specified colormap with specified N on the specified handle
  handle = varargin{1};
  name   = varargin{2};
  n      = varargin{3};

else
  ft_error('incorrect input arguments');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the colormap

% the ones from brewermap are not m-files on disk
brewerlist = repmat(brewermap('list'), [2 1]);
% also include the reverse ones
for k = 1:numel(brewerlist)/2
  brewerlist{k} = sprintf('*%s', brewerlist{k});
end

% the ones from cmocean are not m-files on disk
cmoceanlist = {
  'thermal'
  'haline'
  'solar'
  'ice'
  'gray'
  'oxy'
  'deep'
  'dense'
  'algae'
  'matter'
  'turbid'
  'speed'
  'amp'
  'tempo'
  'rain'
  'phase'
  'topo'
  'balance'
  'delta'
  'curl'
  'diff'
  'tarn'
  };
% also include the reverse ones
for k = 1:numel(cmoceanlist)
  cmoceanlist{end+1} = sprintf('-%s', cmoceanlist{k});
end

% the ones from colorcet are not m-files on disk
colorcetlist = {
  'blueternary'
  'coolwarm'
  'cyclicgrey'
  'depth'
  'divbjy'
  'fire'
  'geographic'
  'geographic2'
  'gouldian'
  'gray'
  'greenternary'
  'grey'
  'heat'
  'phase2'
  'phase4'
  'rainbow'
  'rainbow2'
  'rainbow3'
  'rainbow4'
  'redternary'
  'reducedgrey'
  'yellowheat'
  };
% also include the reverse ones
for k = 1:numel(colorcetlist)
  colorcetlist{end+1} = sprintf('-%s', colorcetlist{k});
end

if isnumeric(name)
  % the user specified an Nx3 array as colormap
  cmap = name;
elseif ismember(name, brewerlist)
  cmap = brewermap(n, name);
elseif startsWith(name, '-') && ismember(name(2:end), brewerlist)
  cmap = brewermap(n, strrep(name, '-', '*'));
elseif ismember(name, cmoceanlist)
  cmap = cmocean(name, n);
elseif ismember(name, colorcetlist)
  cmap = colorcet(name, 'n', n, 'reverse', false);
elseif startsWith(name, '-') && ismember(name(2:end), colorcetlist)
  cmap = colorcet(name(2:end), 'n', n, 'reverse', true);
elseif isequal(name, 'default')
  % requires separate handling, because there's no function called default,
  % the default is taken care of by colormap
  cmap = name;
else
  % this works both for the MATLAB and the MATPLOTLIB colormaps
  % which have the different colormaps available as an m-file
  if name(1)=='-'
    cmap = flipud(feval(name(2:end), n));
  else
    cmap = feval(name, n);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the colormap and/or return it

if nargout>0
  % apply it to the specified figure or axis and return it as a Nx3 array
  colormap(handle, cmap)
else
  % only apply it to the current figure
  colormap(handle, cmap)
  clear cmap
end
