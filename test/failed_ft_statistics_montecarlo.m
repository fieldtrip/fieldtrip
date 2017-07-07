function failed_ft_statistics_montecarlo

% MEM 4gb
% WALLTIME 00:20:00

% TEST test_ft_statistics_montecarlo
% TEST ft_statistics_montecarlo ft_timelockstatistics ft_freqstatistics ft_sourcestatistics clusterstat findcluster

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% test the functionality of ft_statistics_montecarlo, in particular with respect to the clustering behaviour.

% start with some data
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat');
load(filename);

cfg.channel = ft_channelselection('MEG',data.label);
data  = ft_selectdata(cfg, data);

grad  = ft_convert_units(data.grad,'m');
vol   = ft_datatype_headmodel(struct('o',[0 0 0.04],'r',0.08,'unit','m'));

cfg = [];
cfg.dip.pos = [0 sqrt(2)/2 sqrt(2)/2;0 -sqrt(2)/2 sqrt(2)/2]*0.06;
cfg.dip.mom = [1 1;0 0; 0 0];
cfg.grad    = grad;
cfg.vol     = vol;
cfg.channel = 'MEG';
cfg.ntrials = 40;
cfg.absnoise = 8e-8;
cfg.fsample  = 1000;
tmpx = randn(1,1000);
for k = 1:40
  % create some data with an oscillation
  tmp   = ft_preproc_bandpassfilter(randn(1,1000),1000, [15 25]);
  taper = [zeros(1,500) hanning(500)'].*0.8;
  
  tmp2   = ft_preproc_bandpassfilter(tmpx, 1000, [8 12]); 
  taper2 = [hanning(500)',zeros(1,500)];
  
  signal{k} = [tmp.*taper;tmp2.*taper2];
  signal2{k} = zeros(2,1000);
end
cfg.dip.signal = signal;
data1 = ft_dipolesimulation(cfg);
cfg.dip.signal = signal2;
data2 = ft_dipolesimulation(cfg);

cfg = [];
cfg.keeptrials = 'yes';
tlck1 = ft_timelockanalysis(cfg, data1);
tlck2 = ft_timelockanalysis(cfg, data2);

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.design = [ones(1,40) ones(1,40)*2];
cfg.ivar   = 1;
cfg.numrandomization = 50;
stat = ft_timelockstatistics(cfg, tlck1, tlck2);

load ctf275_neighb;
cfg.neighbours = neighbours;
cfg.correctm   = 'cluster';
cfg.clustercritval = [-5 5];
stat = ft_timelockstatistics(cfg, tlck1, tlck2);


%% freq data
cfg = [];
cfg.method = 'mtmconvol';
cfg.keeptrials = 'yes';
cfg.output = 'pow';
cfg.toi    = (0.2:0.05:0.8);
cfg.foi    = 2.5:2.5:40;
cfg.t_ftimwin = ones(1,numel(cfg.foi))*0.4;
cfg.taper  = 'hanning';
freq1      = ft_freqanalysis(cfg, data1);
freq2      = ft_freqanalysis(cfg, data2);

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.design = [ones(1,40) ones(1,40)*2];
cfg.ivar   = 1;
cfg.numrandomization = 50;
cfg.channel = 'MEG';
stat = ft_freqstatistics(cfg, freq1, freq2);

load ctf275_neighb;
cfg.neighbours = neighbours;
cfg.correctm   = 'cluster';
stat = ft_freqstatistics(cfg, freq1, freq2);

%% source data
cfg = [];
cfg.method = 'mtmconvol';
cfg.keeptrials = 'yes';
cfg.output = 'fourier';
cfg.toi    = (0.2:0.05:0.8);
cfg.foi    = 2.5:2.5:40;
cfg.t_ftimwin = ones(1,numel(cfg.foi))*0.4;
cfg.taper  = 'hanning';
freq1      = ft_freqanalysis(cfg, data1);
freq2      = ft_freqanalysis(cfg, data2);
freq       = ft_freqanalysis(cfg, ft_appenddata([], data1, data2));

cfg      = [];
cfg.vol  = vol;
cfg.grad = grad;
cfg.grid.resolution = 0.01;
cfg.channel = 'MEG';
sourcemodel_grid = ft_prepare_leadfield(cfg);

[pnt,tri] = icosahedron642;
cfg.grid  = struct('pos',pnt.*0.06,'tri',tri,'unit','m');
sourcemodel_mesh = ft_prepare_leadfield(cfg);

cfg       = [];
cfg.vol   = vol;
cfg.grid   = sourcemodel_grid;
cfg.method = 'dics';
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda = '10%';
cfg.dics.keepfilter = 'yes';
cfg.dics.realfilter = 'yes';
cfg.frequency = 10;
cfg.latency   = 0.25;
source_grid = ft_sourceanalysis(cfg, freq);

cfg.grid.filter = source_grid.avg.filter;
source1avg = ft_sourceanalysis(cfg, freq1);
source2avg = ft_sourceanalysis(cfg, freq2);

cfg.rawtrial    = 'yes';
source1         = ft_sourceanalysis(cfg, freq1);
source2         = ft_sourceanalysis(cfg, freq2);

cfg.rawtrial    = 'no';
cfg.grid        = sourcemodel_mesh;
source_mesh     = ft_sourceanalysis(cfg, freq);

cfg.grid.filter = source_mesh.avg.filter;
source1avg_mesh = ft_sourceanalysis(cfg, freq1);
source2avg_mesh = ft_sourceanalysis(cfg, freq2);

cfg.rawtrial    = 'yes';
source1_mesh    = ft_sourceanalysis(cfg, freq1);
source2_mesh    = ft_sourceanalysis(cfg, freq2);

%% stats on source data
cfg           = [];
cfg.method    = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.design    = [ones(1,40) ones(1,40)*2];
cfg.ivar      = 1;
cfg.numrandomization = 50;
cfg.parameter = 'pow';
stat          = ft_sourcestatistics(cfg, source1, source2);

% convert the trials into structures with averages, mimicking single
% subject averages
s1 = repmat({source1avg}, [1 40]);
s2 = repmat({source2avg}, [1 40]);
for k = 1:numel(source1.trial)
  s1{k}.avg.pow = source1.trial(k).pow;
  s2{k}.avg.pow = source2.trial(k).pow;
end
cfg.parameter = 'avg.pow';
stat = ft_sourcestatistics(cfg, s1{:}, s2{:});

cfg.correctm = 'cluster';
cfg.parameter = 'pow';
stat = ft_sourcestatistics(cfg, source1, source2);
cfg.parameter = 'avg.pow';
stat = ft_sourcestatistics(cfg, s1{:},   s2{:});

cfg.parameter = 'pow';
cfg.correctm = 'no';
stat = ft_sourcestatistics(cfg, source1_mesh, source2_mesh);

% convert the trials into structures with averages, mimicking single
% subject averages
s1_mesh = repmat({source1avg_mesh}, [1 40]);
s2_mesh = repmat({source2avg_mesh}, [1 40]);
for k = 1:numel(source1_mesh.trial)
  s1_mesh{k}.avg.pow = source1_mesh.trial(k).pow;
  s2_mesh{k}.avg.pow = source2_mesh.trial(k).pow;
end
cfg.parameter = 'avg.pow';
stat = ft_sourcestatistics(cfg, s1_mesh{:}, s2_mesh{:});


cfg.correctm = 'cluster';
cfg.tri      = tri;
cfg.parameter = 'pow';
stat  = ft_sourcestatistics(cfg, source1_mesh, source2_mesh);
cfg.parameter = 'avg.pow';
stat = ft_sourcestatistics(cfg, s1_mesh{:}, s2_mesh{:});


