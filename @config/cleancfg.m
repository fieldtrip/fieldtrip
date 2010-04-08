function [newcfg] = cleancfg(cfg, presentused, defaultused, presentunused, defaultunused);

% CLEANCFG Returns a structure with the config fields that were used
% and displays on screen which fields were used or not.

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

