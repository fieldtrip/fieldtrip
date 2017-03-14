function test_bug2197

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_selectdata

freq = [];
freq.labelcmb = {
  '1' '1'
  '1' '2'
  '1' '3'
  '1' '4'
  '1' '5'
  '1' '6'
  %  '2' '1'
  %  '2' '2'
  '2' '3'
  '2' '4'
  '2' '5'
  '2' '6'
  };
freq.cumtapcnt = ones(15,1);
freq.freq      = 1:7;
freq.crsspctrm = randn(15,10,7);
% freq.crsspctrmdimord = 'rpt_chancmb_freq'; % this is not needed

cfg = [];
cfg.avgoverrpt = 'yes';
output = ft_selectdata(cfg, freq);
assert(size(output.crsspctrm,1)==10)

cfg = [];
cfg.avgoverchancmb = 'yes';
output = ft_selectdata(cfg, freq);
assert(size(output.crsspctrm,2)==1)


cfg = [];
cfg.avgoverrpt = 'yes';
cfg.avgoverchancmb = 'yes';
output = ft_selectdata(cfg, freq);
assert(size(output.crsspctrm,1)==1)
assert(size(output.crsspctrm,2)==7) % frequencies

cfg = [];
cfg.avgoverrpt = 'yes';
cfg.avgoverchancmb = 'yes';
cfg.foilim = [3 5];
output = ft_selectdata(cfg, freq);
assert(size(output.crsspctrm,1)==1)
assert(size(output.crsspctrm,2)==3) % frequencies

