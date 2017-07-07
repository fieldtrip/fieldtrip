function test_bug2727

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_clusterplot

label = {
  'MLC11'    'MLC12'    'MLC13'    'MLC14'    'MLC15'    'MLC21'    'MLC22'    'MLC23'    'MLC24'    'MLC31'    'MLC32'    'MLC33'    'MLC41'    'MLC42' ...
  'MLC43'    'MLF11'    'MLF12'    'MLF21'    'MLF22'    'MLF23'    'MLF31'    'MLF32'    'MLF33'    'MLF34'    'MLF41'    'MLF42'    'MLF43'    'MLF44' ...
  'MLF45'    'MLF51'    'MLF52'    'MLO11'    'MLO12'    'MLO21'    'MLO22'    'MLO31'    'MLO32'    'MLO33'    'MLO41'    'MLO42'    'MLO43'    'MLP11' ...
  'MLP12'    'MLP13'    'MLP21'    'MLP22'    'MLP31'    'MLP32'    'MLP33'    'MLP34'    'MLT11'    'MLT12'    'MLT13'    'MLT14'    'MLT15'    'MLT16' ...
  'MLT21'    'MLT22'    'MLT23'    'MLT24'    'MLT25'    'MLT26'    'MLT31'    'MLT32'    'MLT33'    'MLT34'    'MLT35'    'MLT41'    'MLT42'    'MLT43' ...
  'MLT44'    'MRC11'    'MRC12'    'MRC13'    'MRC14'    'MRC15'    'MRC21'    'MRC22'    'MRC23'    'MRC24'    'MRC31'    'MRC32'    'MRC33'    'MRC41' ...
  'MRC42'    'MRC43'    'MRF11'    'MRF12'    'MRF21'    'MRF22'    'MRF23'    'MRF31'    'MRF32'    'MRF33'    'MRF34'    'MRF41'    'MRF42'    'MRF43' ...
  'MRF44'    'MRF45'    'MRF51'    'MRF52'    'MRO11'    'MRO12'    'MRO21'    'MRO22'    'MRO31'    'MRO32'    'MRO33'    'MRO41'    'MRO42'    'MRO43' ...
  'MRP11'    'MRP12'    'MRP13'    'MRP21'    'MRP22'    'MRP31'    'MRP32'    'MRP33'    'MRP34'    'MRT11'    'MRT12'    'MRT13'    'MRT14'    'MRT15' ...
  'MRT16'    'MRT21'    'MRT22'    'MRT23'    'MRT24'    'MRT25'    'MRT26'    'MRT31'    'MRT32'    'MRT33'    'MRT34'    'MRT35'    'MRT41'    'MRT42' ...
  'MRT43'    'MRT44'    'MZC01'    'MZC02'    'MZF01'    'MZF02'    'MZF03'    'MZO01'    'MZO02'    'MZP01'    'MZP02'
  };

stat = [];
stat.label = label;
stat.prob = rand(151,1);
stat.stat = stat.prob;
stat.dimord = 'chan';
stat.posclusterslabelmat = stat.prob<0.05;
stat.posclusters(1).prob = 0;
% stat.posclusters(1).clusterstat = nan;
% stat.posclusters(1).stddev      = nan;
% stat.posclusters(1).cirange     = nan;

statT = [];
statT.label = label;
statT.prob = rand(151,10);
statT.stat = statT.prob;
statT.time = (1:10)/10; % in seconds
statT.dimord = 'chan_time';
statT.posclusterslabelmat = statT.prob<0.05;
statT.posclusters(1).prob = 0;

statF = [];
statF.label = label;
statF.prob = rand(151,20);
statF.stat = statF.prob;
statF.freq = 1:20; % in Hz
statF.dimord = 'chan_freq';
statF.posclusterslabelmat = statF.prob<0.05;
statF.posclusters(1).prob = 0;

statTF = [];
statTF.label = label;
statTF.prob = rand(151,20,10);
statTF.stat = statTF.prob;
statTF.freq = 1:20; % in Hz
statTF.time = (1:10)/10; % in seconds
statTF.dimord = 'chan_freq_time';
statTF.posclusterslabelmat = statTF.prob<0.05;
statTF.posclusters(1).prob = 0;

statT1F = [];
statT1F.label = label;
statT1F.prob = rand(151,20,1);
statT1F.stat = statT1F.prob;
statT1F.freq = 1:20; % in Hz
statT1F.time = 0.5; % in seconds
statT1F.dimord = 'chan_freq_time';
statT1F.posclusterslabelmat = statT1F.prob<0.05;
statT1F.posclusters(1).prob = 0;

statTF1 = [];
statTF1.label = label;
statTF1.prob = rand(151,1,10);
statTF1.stat = statTF1.prob;
statTF1.freq = 1; % in Hz
statTF1.time = (1:10)/10; % in seconds
statTF1.dimord = 'chan_freq_time';
statTF1.posclusterslabelmat = statTF1.prob<0.05;
statTF1.posclusters(1).prob = 0;

cfg = [];
cfg.layout = 'CTF151.lay';
ft_clusterplot(cfg, stat);
ft_clusterplot(cfg, statT);
ft_clusterplot(cfg, statF);
ft_clusterplot(cfg, statTF1); % single frequency
ft_clusterplot(cfg, statT1F); % single latency

try
  figure; ft_clusterplot(cfg, statTF);
  failed = false;
catch
  failed = true;
end
assert(failed==true, 'this should have failed');




