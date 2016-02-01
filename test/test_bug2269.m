function test_bug2269

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug2296
% TEST ft_preprocessing ft_redefinetrial ft_fetch_data

fs    = 1000;
nsmp  = 1*fs;
nchan = 32;


data = [];
for i=1:nchan
  data.label{i} = sprintf('%d', i);
end
data.time = {[1:nsmp]/fs};
data.trial = {randn(nchan, nsmp)};


begtrl = 1:100:nsmp;
endtrl = begtrl + 99;

cfg = [];
cfg.trl         = [begtrl(:) endtrl(:)];
cfg.trl(:,3)    = 0;
data1 = ft_redefinetrial(cfg, data);


cfg = [];
cfg.trl         = [begtrl(:) endtrl(:)+10]; % now with some overlap
cfg.trl(:,3)    = 0;
cfg.trl(end,:)  = []; % as it goes beyond the end of the data
data2 = ft_redefinetrial(cfg, data);

cfg = [];
cfg.trl = [1 nsmp/2 0];
data1a = ft_redefinetrial(cfg, data1);
data2a = ft_redefinetrial(cfg, data2); % this should not fail, the values are exaclty the same

cfg = [];
cfg.demean = 'yes';
data1b = ft_preprocessing(cfg, data1);
data2b = ft_preprocessing(cfg, data2);

cfg = [];
cfg.trl = [1 nsmp/2 0];
data1c = ft_redefinetrial(cfg, data1b);
try
  data2c = ft_redefinetrial(cfg, data2b); % this should fail, the values are different
  failed = false;
catch
  failed = true;
end

assert(failed, 'ft_redefinetrial should have failed');
