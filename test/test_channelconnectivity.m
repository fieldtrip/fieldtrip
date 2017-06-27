function test_channelconnectivity

% MEM 1000mb
% WALLTIME 00:10:00

% TEST test_channelconnectivity channelconnectivity ft_megplanar statistics_montecarlo

% Tests the functionality of private/channelconnectivity(), which generates
% a logical NxN matrix of channel connectivity based on a neighbours
% structure.

fprintf('testing the following version of channelconnectivity():\n%s\n', which('channelconnectivity'));

%% generate random 'true' connectivity and channel labels
nchan = 100;
conn = false(nchan,nchan);
labels = {};

symbols = ['a':'z' 'A':'Z' '0':'9'];
for k = 1:nchan
  % every channel has between 0 and 9 random neighbours
  conn(k,randperm(nchan,randi(10))) = 1;
  
  % random channel name length between 1 and 7
  labels{end+1} = symbols(randi(numel(symbols),1,randi(7)));
  while numel(unique(labels)) ~= numel(labels) % ensure unique channel names
    labels{end} = symbols(randi(numel(symbols),1,randi(7)));
  end
end

%% generate 'neighbours' structure
neighb = [];
for k = 1:nchan
  neighb(k).label = labels{k};
  neighb(k).neighblabel = labels(conn(k,:))';
end

%% compute connectivity from neighbours and verify

% case 1: channels in cfg, all channels present in cfg and neighbours
fprintf('test case 1:');
cfg = [];
cfg.channel = labels;
cfg.neighbours = neighb;
connA = channelconnectivity(cfg);
if ~isequal(conn, connA)
  error('estimated neighbourhood connectivity wrong at case 1');
end
fprintf('success\n');

% case 2: channels in cfg, some channels missing from cfg.channel
fprintf('test case 2:');
cfg = [];
inds = randperm(numel(labels), randi(nchan));
cfg.channel = labels(inds);
cfg.neighbours = neighb;
connA = channelconnectivity(cfg);
if ~isequal(conn(inds,inds), connA)
  error('estimated neighbourhood connectivity wrong at case 2');
end
fprintf('success\n');

% case 3: channels in cfg, some channels missing from cfg.neighbours
fprintf('test case 3:');
cfg = [];
inds = randperm(numel(labels), randi([nchan/2 nchan]));
cfg.channel = labels;
cfg.neighbours = neighb(inds);
connA = channelconnectivity(cfg);
allinds = true(nchan,1);
allinds(inds) = 0;
missingconn = connA(allinds,:); % these should all be 0 because neighbour information was not there for this channel
if ~isequal(conn(inds,:), connA(inds,:)) || any(missingconn(:))
  error('estimated neighbourhood connectivity wrong at case 3');
end
fprintf('success\n');

% case 4: channels in cfg, some channels missing from both cfg.channel and cfg.neighbours
fprintf('test case 4:');
cfg = [];
inds1 = randperm(numel(labels), randi(nchan));
inds2 = randperm(numel(labels), randi([nchan/2 nchan]));
cfg.channel = labels(inds1);
cfg.neighbours = neighb(inds2);
connA = channelconnectivity(cfg);
% FIXME: not sure how to check whether the answer here is correct
fprintf('success\n');

% case 5: channels in data, all channels present in cfg and neighbours
fprintf('test case 5:');
cfg = [];
data = [];
data.label = labels;
cfg.neighbours = neighb;
connA = channelconnectivity(cfg, data);
if ~isequal(conn, connA)
  error('estimated neighbourhood connectivity wrong at case 5');
end
fprintf('success\n');

% case 6: channels in data, some channels missing from cfg.channel
fprintf('test case 6:');
cfg = [];
inds = randperm(numel(labels), randi(nchan));
data = [];
data.label = labels(inds);
cfg.neighbours = neighb;
connA = channelconnectivity(cfg, data);
if ~isequal(conn(inds,inds), connA)
  error('estimated neighbourhood connectivity wrong at case 6');
end
fprintf('success\n');

% case 7: channels in data, some channels missing from cfg.neighbours
fprintf('test case 7:');
cfg = [];
inds = randperm(numel(labels), randi([nchan/2 nchan]));
data = [];
data.label = labels;
cfg.neighbours = neighb(inds);
connA = channelconnectivity(cfg, data);
allinds = true(nchan,1);
allinds(inds) = 0;
missingconn = connA(allinds,:); % these should all be 0 because neighbour information was not there for this channel
if ~isequal(conn(inds,:), connA(inds,:)) || any(missingconn(:))
  error('estimated neighbourhood connectivity wrong at case 7');
end
fprintf('success\n');

% case 8: channels in data, some channels missing from both cfg.channel and cfg.neighbours
fprintf('test case 8:');
cfg = [];
inds1 = randperm(numel(labels), randi(nchan));
inds2 = randperm(numel(labels), randi([nchan/2 nchan]));
data = [];
data.label = labels(inds1);
cfg.neighbours = neighb(inds2);
connA = channelconnectivity(cfg, data);
% FIXME: not sure how to check whether the answer here is correct
fprintf('success\n');

end
