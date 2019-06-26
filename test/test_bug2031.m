function test_bug2031

% MEM 1500mb
% WALLTIME 00:10:00

% DEPENDENCY ft_topoplotER ft_topoplotTFR ft_singleplotER ft_singleplotTFR

% this test script makes some figures, to verify that the FieldTrip menu
% item is not added multiple times (1) when subplots are being used, and
% (2) when plotting in a figure that already exists.

%% make some data
tl = [];
tl.label = {'MLC11','MLC12','MLC13','MLC14','MLC15','MLC16','MLC17','MLC21','MLC22','MLC23','MLC24','MLC25','MLC31','MLC32','MLC41','MLC42','MLC51','MLC52','MLC53','MLC54','MLC55','MLC61','MLC62','MLC63','MLF11','MLF12','MLF13','MLF14','MLF21','MLF22','MLF23','MLF24','MLF25','MLF31','MLF32','MLF33','MLF34','MLF35','MLF41','MLF42','MLF43','MLF44','MLF45','MLF46','MLF51','MLF52','MLF53','MLF54','MLF55','MLF56','MLF61','MLF62','MLF63','MLF64','MLF65','MLF66','MLF67','MLO13','MLO14','MLO23','MLO24','MLO31','MLO32','MLO33','MLO34','MLO41','MLO42','MLO43','MLO44','MLO51','MLO52','MLO53','MLP11','MLP12','MLP21','MLP22','MLP23','MLP31','MLP32','MLP33','MLP34','MLP35','MLP41','MLP42','MLP43','MLP44','MLP45','MLP51','MLP52','MLP53','MLP54','MLP55','MLP56','MLP57','MLT11','MLT12','MLT13','MLT14','MLT15','MLT16','MLT21','MLT22','MLT23','MLT24','MLT25','MLT26','MLT27','MLT31','MLT32','MLT33','MLT34','MLT35','MLT36','MLT37','MLT41','MLT42','MLT43','MLT44','MLT45','MLT46','MLT47','MLT51','MLT52','MLT53','MLT54','MLT55','MLT56','MLT57','MRC11','MRC12','MRC13','MRC14','MRC15','MRC16','MRC17','MRC21','MRC22','MRC23','MRC24','MRC25','MRC31','MRC32','MRC41','MRC42','MRC51','MRC52','MRC53','MRC54','MRC55','MRC61','MRC62','MRC63','MRF11','MRF12','MRF13','MRF14','MRF21','MRF22','MRF23','MRF24','MRF25','MRF31','MRF32','MRF33','MRF34','MRF35','MRF41','MRF42','MRF43','MRF44','MRF45','MRF46','MRF51','MRF52','MRF53','MRF54','MRF55','MRF56','MRF61','MRF62','MRF63','MRF64','MRF65','MRF67','MRO11','MRO12','MRO13','MRO14','MRO21','MRO22','MRO23','MRO24','MRO31','MRO32','MRO33','MRO34','MRO41','MRO42','MRO43','MRO44','MRO51','MRO53','MRP11','MRP12','MRP21','MRP22','MRP23','MRP31','MRP32','MRP33','MRP34','MRP35','MRP41','MRP42','MRP43','MRP44','MRP45','MRP51','MRP52','MRP53','MRP54','MRP55','MRP56','MRP57','MRT11','MRT12','MRT13','MRT14','MRT15','MRT16','MRT21','MRT22','MRT23','MRT24','MRT25','MRT26','MRT27','MRT31','MRT32','MRT33','MRT34','MRT35','MRT36','MRT37','MRT41','MRT42','MRT43','MRT44','MRT45','MRT46','MRT47','MRT51','MRT52','MRT53','MRT54','MRT55','MRT56','MRT57','MZC01','MZC02','MZC03','MZC04','MZF01','MZF02','MZF03','MZO01','MZO03','MZP01','MLO11','MLO12','MLO21','MLO22','MZO02'};
tl.time = 0.01:0.01:1;
tl.dimord = 'chan_time';
tl.avg = randn(numel(tl.label), numel(tl.time));

freq = [];
freq.label = tl.label;
freq.freq = 1:30;
freq.time = 0.01:0.01:1;
freq.dimord = 'chan_freq_time';
freq.powspctrm = randn(numel(freq.label), numel(freq.freq), numel(freq.time));
  
%% ft_topoplotER, subplots
% expected behavior: no menu after second and subsequent calls

figure;
cfg = [];
cfg.layout = 'CTF275.lay';
subplot(2,2,1);
ft_topoplotER(cfg, tl);
subplot(2,2,2);
ft_topoplotER(cfg, tl);
subplot(2,2,3);
ft_topoplotER(cfg, tl);
subplot(2,2,4);
ft_topoplotER(cfg, tl);

%% ft_topoplotTFR, subplots
% expected behavior: no menu after second and subsequent calls

figure;
cfg = [];
cfg.layout = 'CTF275.lay';
subplot(2,2,1);
ft_topoplotTFR(cfg, freq);
subplot(2,2,2);
ft_topoplotTFR(cfg, freq);
subplot(2,2,3);
ft_topoplotTFR(cfg, freq);
subplot(2,2,4);
ft_topoplotTFR(cfg, freq);

%% ft_topoplotER, subplots
% expected behavior: no menu after second and subsequent calls

figure;
cfg = [];
cfg.layout = 'CTF275.lay';
subplot(2,2,1);
ft_topoplotER(cfg, tl);
subplot(2,2,2);
ft_topoplotER(cfg, tl);
subplot(2,2,3);
ft_topoplotER(cfg, tl);
subplot(2,2,4);
ft_topoplotER(cfg, tl);

%% ft_topoplotER, subsequent plots
% expected behavior: one menu after second and subsequent calls

figure;
cfg = [];
cfg.layout = 'CTF275.lay';
ft_topoplotER(cfg, tl);
ft_topoplotER(cfg, tl);
ft_topoplotER(cfg, tl);

%% ft_topoplotTFR, subsequent plots
% expected behavior: one menu after second and subsequent calls

figure;
cfg = [];
cfg.layout = 'CTF275.lay';
ft_topoplotTFR(cfg, freq);
ft_topoplotTFR(cfg, freq);
ft_topoplotTFR(cfg, freq);

%% ft_singleplotER, subplots
% expected behavior: no menu after second and subsequent calls

figure;
cfg = [];
subplot(2,2,1);
ft_singleplotER(cfg, tl);
subplot(2,2,2);
ft_singleplotER(cfg, tl);
subplot(2,2,3);
ft_singleplotER(cfg, tl);
subplot(2,2,4);
ft_singleplotER(cfg, tl);

%% ft_singleplotTFR, subplots
% expected behavior: no menu after second and subsequent calls

figure;
cfg = [];
subplot(2,2,1);
ft_singleplotTFR(cfg, freq);
subplot(2,2,2);
ft_singleplotTFR(cfg, freq);
subplot(2,2,3);
ft_singleplotTFR(cfg, freq);
subplot(2,2,4);
ft_singleplotTFR(cfg, freq);

%% ft_singleplotER, subsequent plots
% expected behavior: one menu after second and subsequent calls

figure;
cfg = [];
ft_singleplotER(cfg, tl);
ft_singleplotER(cfg, tl);
ft_singleplotER(cfg, tl);

%% ft_singleplotTFR, subsequent plots
% expected behavior: one menu after second and subsequent calls

figure;
cfg = [];
ft_singleplotTFR(cfg, freq);
ft_singleplotTFR(cfg, freq);
ft_singleplotTFR(cfg, freq);
