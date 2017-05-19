function test_bug3295

% WALLTIME 00:10:00
% MEM 1gb

n_chans = 10;

%%
% create a raw data structure, 2 channels each....

data_2 = {};
for i = 1:n_chans/2
  cfg = [];
  cfg.method = 'broadband';
  cfg.fsample = 100;
  cfg.trllen = 2;
  cfg.numtrl = 50;
  
  cfg.n1.ampl = 1;
  cfg.n1.bpfreq = [5 10];
  cfg.n2.ampl = 2;
  cfg.n2.bpfreq = [12 20];
  cfg.noise.ampl = 3;
  
  tmp_data = ft_freqsimulation(cfg);
  
  cfg = [];
  cfg.channel = {'n1', 'n2'};
  
  tmp_data = ft_selectdata(cfg, tmp_data);
  tmp_data.label{1} = sprintf('chan_%02d', (i-1)*2);
  tmp_data.label{2} = sprintf('chan_%02d', (i-1)*2 + 1);
  
  data_2{end+1} = tmp_data;
end %for

%%
% append raw data, they all have different channels

cfg = [];
cfg.appenddim = 'chan';
data_all = ft_appenddata(cfg, data_2{:});

assert(numel(data_all.label)==10);

cfg = [];
% cfg.appenddim = 'chan'; % should be auto-detected
data_all = ft_appenddata(cfg, data_2{:});

assert(numel(data_all.label)==10);

%%
% create two other dataset with fewer channels...

cfg = [];
cfg.channel = 1:7;

data_few_1 = ft_selectdata(cfg, data_all);

cfg = [];
cfg.channel = 5:10;

data_few_2 = ft_selectdata(cfg, data_all);

%%
% do appenddata, this should return the intersection of channels

cfg = [];
cfg.appenddim = 'rpt'; % auto-detection does not work, since some - but not all - channels overlap
data_final = ft_appenddata(cfg, data_all, data_few_1, data_few_2);

assert(numel(data_final.label)==3);

%%
% do selectdata to get correct intersection of channels...

cfg = [];
[data_all_selected, data_few_1_selected, data_few_2_selected] = ft_selectdata(cfg, data_all, data_few_1, data_few_2);

assert(numel(data_all_selected.label)==3);
assert(numel(data_few_1_selected.label)==3);
assert(numel(data_few_2_selected.label)==3);

cfg = [];
cfg.appenddim = []; % auto-detection should work
data_compare = ft_appenddata(cfg, data_all_selected, data_few_1_selected, data_few_2_selected);

assert(isequal(keepfields(data_compare, {'trial', 'time', 'label'}), keepfields(data_final, {'trial', 'time', 'label'})))
