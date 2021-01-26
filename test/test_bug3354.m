function test_bug3354

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_selectdata

trial = cell(1,10);
time  = cell(1,10);
for k = 1:10
  trial{k} = randn(2,1000);
  time{k}  = (1:1000)./1000;
end
data.trial = trial;
data.time  = time;
data.label = {'chan01';'chan02'};

data = ft_datatype_raw(data);
data = ft_checkdata(data, 'hassampleinfo', 'yes');

cfg = [];
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.latency = [0.20 0.699];
data2 = ft_selectdata(cfg, data);
tlck2 = ft_selectdata(cfg, tlck);

assert(isfield(data2, 'sampleinfo'));
assert(isfield(tlck2, 'sampleinfo'));

for k = 1:numel(data2.trial)
  tmp1 = data.time{k}(data2.sampleinfo(k,1)-data.sampleinfo(k,1)+1);
  tmp2 = data2.time{k}(1);
  assert(isequal(tmp1, tmp2));
end

tmp1 = tlck.time(tlck2.sampleinfo(:,1)-tlck.sampleinfo(:,1)+1);
tmp2 = tlck2.time(1);
assert(all(tmp1==tmp2));
