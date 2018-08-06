function [newcfg] = cleancfg(cfg, presentused, defaultused, presentunused, defaultunused)

% CLEANCFG Returns a structure with the config fields that were used
% and displays on screen which fields were used or not.

% Copyright (C) 2012-2015, Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
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

% set the defaults
if nargin<2
  presentused = 0;
end
if nargin<3
  defaultused = 0;
end
if nargin<4
  presentunused = 0;
end
if nargin<5
  defaultunused = 0;
end

r = access(cfg, 'reference');
a = access(cfg, 'assign');
o = access(cfg, 'original');
v = access(cfg, 'value');

key      = fieldnames(cfg); key = key(:)';
used     = zeros(size(key));
modified = zeros(size(key));
original = zeros(size(key));

for i=1:length(key)
  used(i)     = (r.(key{i})>0);
  modified(i) = (a.(key{i})>0);
  original(i) = (o.(key{i})>0);
end

if presentused
  fprintf('\nThe following config fields were USED as specified by you\n');
  sel = find(used & ~modified);
  if numel(sel)
    fprintf('  cfg.%s\n', key{sel});
  else
    fprintf('  <none>\n');
  end
end

if defaultused
  fprintf('\nThe following config fields were USED and were adjusted\n');
  sel = find(used & modified);
  if numel(sel)
    fprintf('  cfg.%s\n', key{sel});
  else
    fprintf('  <none>\n');
  end
end

if presentunused
  fprintf('\nThe following config fields were NOT USED and were specified by you\n');
  sel = find(~used & original);
  if numel(sel)
    fprintf('  cfg.%s\n', key{sel});
  else
    fprintf('  <none>\n');
  end
end

if defaultunused
  fprintf('\nThe following config fields were NOT USED and set to defaults\n');
  sel = find(~used & ~original);
  if numel(sel)
    fprintf('  cfg.%s\n', key{sel});
  else
    fprintf('  <none>\n');
  end
end

if nargout
  usedkey = key(find(used));
  usedval = {};
  for i=1:length(usedkey)
    usedval{i} = v.(usedkey{i});
  end
  newcfg = cell2struct(usedval, usedkey, 2);
end

