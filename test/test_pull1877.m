function test_pull1877

% MEM 1gb
% WALLTIME 00:20:00
% DEPENDENCY ft_prepare_neighbours ft_sourceplot ft_checkdata

% make a parcellation:
[pos, tri] = mesh_sphere(1026);

atlas.pos = pos;
atlas.tri = tri;
atlas.parcellation = zeros(size(atlas.pos,1),1);
atlas.parcellation(atlas.pos(:,1)<=0 & atlas.pos(:,2)<=0 & atlas.pos(:,3)<=0) = 1;
atlas.parcellation(atlas.pos(:,1)>0  & atlas.pos(:,2)<=0 & atlas.pos(:,3)<=0) = 4;
atlas.parcellation(atlas.pos(:,1)<=0 & atlas.pos(:,2)>0  & atlas.pos(:,3)<=0) = 2;
atlas.parcellation(atlas.pos(:,1)>0  & atlas.pos(:,2)>0  & atlas.pos(:,3)<=0) = 3;
atlas.parcellation(atlas.pos(:,1)<=0 & atlas.pos(:,2)<=0 & atlas.pos(:,3)>0) = 5;
atlas.parcellation(atlas.pos(:,1)>0  & atlas.pos(:,2)<=0 & atlas.pos(:,3)>0) = 8;
atlas.parcellation(atlas.pos(:,1)<=0 & atlas.pos(:,2)>0  & atlas.pos(:,3)>0) = 6;
atlas.parcellation(atlas.pos(:,1)>0  & atlas.pos(:,2)>0  & atlas.pos(:,3)>0) = 7;
atlas.parcellationlabel ={'1';'2';'3';'4';'5';'6';'7';'8'};

cfg = [];
cfg.method = 'parcellation';
neighbours = ft_prepare_neighbours(cfg, atlas);

% this is expected per the above specification
assert(isequal(neighbours(1).neighblabel, {'2' '4' '5'}));
assert(isequal(neighbours(2).neighblabel, {'1' '3' '6'}));
assert(isequal(neighbours(3).neighblabel, {'2' '4' '7'}));
assert(isequal(neighbours(4).neighblabel, {'1' '3' '8'}));
assert(isequal(neighbours(5).neighblabel, {'1' '6' '8'}));
assert(isequal(neighbours(6).neighblabel, {'2' '5' '7'}));
assert(isequal(neighbours(7).neighblabel, {'3' '6' '8'}));
assert(isequal(neighbours(8).neighblabel, {'4' '5' '7'}));

% make some data
data = [];
data.trial = randn(200, 8, 1);
data.trial(1:100, [1 2], :) = data.trial(1:100, [1 2], :) + 2;
data.trial(1:100, [6 7], :) = data.trial(1:100, [6 7], :) - 2;
data.dimord = 'rpt_chan_time';
data.time = 1;
data.label = atlas.parcellationlabel;

cfg           = [];
cfg.method    = 'montecarlo';
cfg.numrandomization = 1000;
cfg.design    = [ones(1,100) ones(1,100)*2];
cfg.statistic = 'indepsamplesT';
cfg.correctm  = 'cluster';
cfg.neighbours = neighbours;
cfg.clusteralpha = 0.0001;
stat = ft_timelockstatistics(cfg, data);

assert(isequal(stat.posclusterslabelmat, [1 1 0 0 0 0 0 0]'));
assert(isequal(stat.negclusterslabelmat, [0 0 0 0 0 1 1 0]'));

% check whether swapping the labels in the data has an effect on the
% clustering etc:
reorder = [1 3 2 4 5 6 8 7];
data2 = data;
data2.label = data.label(reorder);
data2.trial = data.trial(:, reorder, :);
stat2 = ft_timelockstatistics(cfg, data2);

assert(isequal(stat2.posclusterslabelmat, [1 0 1 0 0 0 0 0]'));
assert(isequal(stat2.negclusterslabelmat, [0 0 0 0 0 1 0 1]'));

stat.brainordinate = atlas;

cfg = [];
cfg.funparameter = 'stat';
cfg.method = 'surface';
ft_sourceplot(cfg, stat);

stat2.brainordinate = atlas;
ft_sourceplot(cfg, stat2);
