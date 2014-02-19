function test_bug2069

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2069 ft_appendfreq ft_freqstatistics ft_timelockstatistics

%% simulate some data

tl = cell(1,20);
tl_shuff = cell(1,20);
freq = cell(1,20);
freq_shuff = cell(1,20);

for k = 1:20
  % tl will hold timelocked structures with no channel order shuffling
  tl{k} = [];
  tl{k}.label = {'MLC11','MLC12','MLC13','MLC14','MLC15','MLC16','MLC17','MLC21','MLC22','MLC23','MLC24','MLC25','MLC31','MLC32','MLC41','MLC42','MLC51','MLC52','MLC53','MLC54','MLC55','MLC61','MLC62','MLC63','MLF11','MLF12','MLF13','MLF14','MLF21','MLF22','MLF23','MLF24','MLF25','MLF31','MLF32','MLF33','MLF34','MLF35','MLF41','MLF42','MLF43','MLF44','MLF45','MLF46','MLF51','MLF52','MLF53','MLF54','MLF55','MLF56','MLF61','MLF62','MLF63','MLF64','MLF65','MLF66','MLF67','MLO13','MLO14','MLO23','MLO24','MLO31','MLO32','MLO33','MLO34','MLO41','MLO42','MLO43','MLO44','MLO51','MLO52','MLO53','MLP11','MLP12','MLP21','MLP22','MLP23','MLP31','MLP32','MLP33','MLP34','MLP35','MLP41','MLP42','MLP43','MLP44','MLP45','MLP51','MLP52','MLP53','MLP54','MLP55','MLP56','MLP57','MLT11','MLT12','MLT13','MLT14','MLT15','MLT16','MLT21','MLT22','MLT23','MLT24','MLT25','MLT26','MLT27','MLT31','MLT32','MLT33','MLT34','MLT35','MLT36','MLT37','MLT41','MLT42','MLT43','MLT44','MLT45','MLT46','MLT47','MLT51','MLT52','MLT53','MLT54','MLT55','MLT56','MLT57','MRC11','MRC12','MRC13','MRC14','MRC15','MRC16','MRC17','MRC21','MRC22','MRC23','MRC24','MRC25','MRC31','MRC32','MRC41','MRC42','MRC51','MRC52','MRC53','MRC54','MRC55','MRC61','MRC62','MRC63','MRF11','MRF12','MRF13','MRF14','MRF21','MRF22','MRF23','MRF24','MRF25','MRF31','MRF32','MRF33','MRF34','MRF35','MRF41','MRF42','MRF43','MRF44','MRF45','MRF46','MRF51','MRF52','MRF53','MRF54','MRF55','MRF56','MRF61','MRF62','MRF63','MRF64','MRF65','MRF67','MRO11','MRO12','MRO13','MRO14','MRO21','MRO22','MRO23','MRO24','MRO31','MRO32','MRO33','MRO34','MRO41','MRO42','MRO43','MRO44','MRO51','MRO53','MRP11','MRP12','MRP21','MRP22','MRP23','MRP31','MRP32','MRP33','MRP34','MRP35','MRP41','MRP42','MRP43','MRP44','MRP45','MRP51','MRP52','MRP53','MRP54','MRP55','MRP56','MRP57','MRT11','MRT12','MRT13','MRT14','MRT15','MRT16','MRT21','MRT22','MRT23','MRT24','MRT25','MRT26','MRT27','MRT31','MRT32','MRT33','MRT34','MRT35','MRT36','MRT37','MRT41','MRT42','MRT43','MRT44','MRT45','MRT46','MRT47','MRT51','MRT52','MRT53','MRT54','MRT55','MRT56','MRT57','MZC01','MZC02','MZC03','MZC04','MZF01','MZF02','MZF03','MZO01','MZO03','MZP01','MLO11','MLO12','MLO21','MLO22','MZO02'};
  
  tl{k}.avg = randn(273,100);
  tl{k}.time = 0.01:0.01:1;
  tl{k}.dimord = 'chan_time';
  
  % shuffle the channel orders for tl_shuff
  perm = randperm(numel(tl{k}.label));
  tl_shuff{k} = tl{k};
  tl_shuff{k}.label = tl{k}.label(perm);
  tl_shuff{k}.avg = tl{k}.avg(perm,:);

  % take the same data as tl for a freq structure with singleton freq
  freq{k} = [];
  freq{k}.label = tl{k}.label;
  freq{k}.freq = 10;
  freq{k}.time = tl{k}.time;
  freq{k}.powspctrm = reshape(tl{k}.avg, [273 1 100]);
  freq{k}.dimord = 'chan_freq_time';
  
  % and also make a freq_shuff using the same shuffling as for tl_shuff
  freq_shuff{k} = freq{k};
  freq_shuff{k}.label = freq{k}.label(perm);
  freq_shuff{k}.powspctrm = freq{k}.powspctrm(perm,:,:);
end

% note: to visually verify that the shuffling was done correctly simply do
% cfg = [];
% cfg.layout = 'CTF275.lay';
% figure; ft_topoplotER(cfg, tl{1});
% figure; ft_topoplotER(cfg, tl_shuff{1});

%% stats cfg: common stuff for both methods
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;

c2 = [];
c2.method = 'template';
c2.template = 'CTF275_neighb.mat';
cfg.neighbours = ft_prepare_neighbours(c2);

cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 10;
cfg.clusterthreshold = 'nonparametric_common';

cfg.design = [
  ones(1,10) ones(1,10)*2;
  1:10 1:10
];
cfg.ivar  = 1;
cfg.uvar = 2;

%% do the stats

% timelockstatistics
stat_tl = ft_timelockstatistics(cfg, tl{:});
stat_tl_shuff = ft_timelockstatistics(cfg, tl_shuff{:});

% freqstatistics
cfg.avgoverfreq = 'yes';
cfg.frequency = [7 14];
stat_freq = ft_freqstatistics(cfg, freq{:});
stat_freq_shuff = ft_freqstatistics(cfg, freq_shuff{:});

%% diagnostics

% make channel order across structures consistent again
% note that we only fix the .stat field as that is the only one we test
% later on
[sel1,sel2] = match_str(stat_tl.label, stat_tl_shuff.label);
stat_tl_shuff.label = stat_tl_shuff.label(sel2);
stat_tl_shuff.stat = stat_tl_shuff.stat(sel2,:);
stat_tl_shuff.prob = stat_tl_shuff.prob(sel2,:);

[sel1,sel2] = match_str(stat_tl.label, stat_freq.label);
stat_freq.label = stat_freq.label(sel2);
stat_freq.stat = stat_freq.stat(sel2,:,:);
stat_freq.prob = stat_freq.prob(sel2,:,:);

[sel1,sel2] = match_str(stat_tl.label, stat_freq_shuff.label);
stat_freq_shuff.label = stat_freq_shuff.label(sel2);
stat_freq_shuff.stat = stat_freq_shuff.stat(sel2,:,:);
stat_freq_shuff.prob = stat_freq_shuff.prob(sel2,:,:);

fprintf('\ndiagnostic results:\n');
errflag = 0;

if all(stat_tl.stat(:) == stat_tl_shuff.stat(:))
  fprintf('for timelockstatistics, shuffling does not matter\n');
else
  fprintf('ERROR in timelockstatistics when channels are shuffled!\n');
  errflag = 1;
end

if all(stat_freq.stat(:) == stat_freq_shuff.stat(:))
  fprintf('for freqstatistics, shuffling does not matter\n');
else
  fprintf('ERROR in freqstatistics when channels are shuffled!\n');
  errflag = 1;
end

if all(stat_freq.stat(:) == stat_tl.stat(:))
  fprintf('with normal channel order, freqstats equal tlstats\n');
else
  fprintf('with normal channel order, freqstats DOES NOT equal tlstats\n');
  errflag = 1;
end

if all(stat_freq_shuff.stat(:) == stat_tl_shuff.stat(:))
  fprintf('with identically shuffled channel order, freqstats equal tlstats\n');
else
  fprintf('with identically shuffled channel order, freqstats DOES NOT equal tlstats\n');
  errflag = 1;
end

if errflag
  fprintf('\n');
  error('something bad happened, see above for details');
end
  
  
  
  
  
