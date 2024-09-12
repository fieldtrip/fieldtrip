function test_ft_baddata

% WALLTIME 00:10:00
% DATA public
% MEM 1gb

%%

[ftver, ftpath] = ft_version;
gradfile = fullfile(ftpath, 'template', 'gradiometer', 'ctf151.mat');

grad = ft_read_sens(gradfile);

%%
% compute simulated data with a realistic amplitude and realistic noise
% later we will add artifacts to this data

chansel = strcmp(grad.chantype, 'meggrad');

grad.tra      = grad.tra(chansel,:);
grad.chanpos  = grad.chanpos(chansel,:);
grad.chanori  = grad.chanori(chansel,:);
grad.label    = grad.label(chansel);
grad.chantype = grad.chantype(chansel);
grad.chanunit = grad.chanunit(chansel);

%%

headmodel = [];
headmodel.type = 'singlesphere';
headmodel.unit = 'm';
headmodel.r = 0.09;
headmodel.o = [0.02 0 0.04];

%%

if false
  figure
  hold on
  ft_plot_headmodel(headmodel, 'facealpha', 0.5, 'facecolor', 'skin');
  ft_plot_sens(grad);
  ft_plot_axes(grad);
end

%%

cfg = [];
cfg.sourcemodel.unit = 'm';
cfg.sourcemodel.pos = [0 0 0.08];
cfg.sourcemodel.mom = [100 0 0] * 1e-9;
cfg.grad = grad;
cfg.headmodel = headmodel;
cfg.relnoise = 0;
cfg.absnoise = 100e-15;
cfg.sourcemodel.frequency = 10;
cfg.sourcemodel.phase = 0;
cfg.sourcemodel.amplitude = 1;
data = ft_dipolesimulation(cfg);

% add sampleinfo
data = ft_checkdata(data, 'hassampleinfo', 'yes');

%%

cfg = [];
cfg.layout = 'CTF151_helmet';
ft_databrowser(cfg, data)

%%

ntrl  = length(data.trial);
nchan = length(data.label);
ntime = length(data.time{1});

% large variance
badchannel = 1; trl = 1;
data.trial{trl}(badchannel,:) = data.trial{trl}(badchannel,:) + randn(1,ntime) * 1e-12;
badchannel = 2; trl = 2;
data.trial{trl}(badchannel,:) = data.trial{trl}(badchannel,:) + randn(1,ntime) * 1e-12;

% large positive offset
badchannel = 3; trl = 3;
data.trial{trl}(badchannel,:) = data.trial{trl}(badchannel,:) + 10 * 1e-12;
badchannel = 4; trl = 4;
data.trial{trl}(badchannel,:) = data.trial{trl}(badchannel,:) + 10 * 1e-12;

% large negative offset
badchannel = 5; trl = 5;
data.trial{trl}(badchannel,:) = data.trial{trl}(badchannel,:) - 10 * 1e-12;
badchannel = 6; trl = 6;
data.trial{trl}(badchannel,:) = data.trial{trl}(badchannel,:) - 10 * 1e-12;

% small variance
badchannel = 7; trl = 7;
data.trial{trl}(badchannel,:) = 0.01 * data.trial{trl}(badchannel,:);
badchannel = 8; trl = 8;
data.trial{trl}(badchannel,:) = 0.01 * data.trial{trl}(badchannel,:);

%%


cfg = [];
cfg.method = 'template';
neighbours = ft_prepare_neighbours(cfg, data);

%%

cfg = [];
cfg.metric = 'std';
cfg.threshold = 0.5 * 1e-12;
tmpcfg = ft_badchannel(cfg, data);
% channel 1 and 2 have large variance
assert(contains(data.label{1}, tmpcfg.badchannel));
assert(contains(data.label{2}, tmpcfg.badchannel));

%%

cfg = [];
cfg.metric = 'std';
cfg.threshold = 0.5 * 1e-12;
tmpcfg = ft_badsegment(cfg, data);
% trial 1 and 2 have large variance
assert(isequal(tmpcfg.artfctdef.badsegment.artifact(1,:), data.sampleinfo(1,:)));
assert(isequal(tmpcfg.artfctdef.badsegment.artifact(1,:), data.sampleinfo(1,:)));

%%

cfg = [];
cfg.metric = 'std';
cfg.threshold = 0.5 * 1e-12;
data_clean = ft_baddata(cfg, data);
% channel 1 in trial 1 has large variance
% channel 2 in trial 2 has large variance
assert(all(isnan(data_clean.trial{1}(1,:))))
assert(all(isnan(data_clean.trial{2}(2,:))))
for i=3:ntrl
  for j=1:nchan
    assert(~any(isnan(data_clean.trial{i}(j,:))))
  end
end

%%

if false
  cfg = [];
  cfg.layout = 'CTF151_helmet';
  ft_databrowser(cfg, data);

  cfg = [];
  cfg.layout = 'CTF151_helmet';
  ft_databrowser(cfg, data_clean);

  cfg = [];
  cfg.layout = 'CTF151_helmet';
  cfg.neighbours = neighbours;
  ft_rejectvisual(cfg, data);

  cfg = [];
  cfg.layout = 'CTF151_helmet';
  cfg.neighbours = neighbours;
  ft_rejectvisual(cfg, data_clean)
end
