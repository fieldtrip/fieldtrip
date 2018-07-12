function test_ft_denoise_tsr

% MEM 8000mb
% WALLTIME 01:30:00

% TEST test_ft_denoise_tsr
% TEST ft_denoise_tsr

% create some data
for k = 1:50
  trial{1,k} = ft_preproc_bandpassfilter(randn(1,1000), 1000, [8 12], [], 'firws');
  time{1,k}  = (0:999)./1000;
end

refdata.trial = trial;
refdata.time  = time;
refdata.label = {'refchan'};

krn = [zeros(1,20) linspace(1,0,20)];
for k = 1:50
  trial{k} = convn(trial{k},krn,'same');
end
data = refdata;
data.trial = trial;
data.label = {'chan01'};

cfg = [];
cfg.refchannel = {'refchan'};
cfg.reflags    = (-10:30)./1000;
cfg.method     = 'mlr';
% out1 = ft_denoise_tsr(cfg, data, refdata);
% 
% cfg.method     = 'mlrridge';
% cfg.threshold  = [0.1 0];
% out2 = ft_denoise_tsr(cfg, data, refdata);

for k = 1:numel(trial)
  trial{k} = trial{k} + 10.*ft_preproc_bandpassfilter(randn(1,1000), 1000, [8 12], [], 'firws');
end
data.trial = trial;

cfg = [];
cfg.refchannel = {'refchan'};
cfg.reflags    = (-10:30)./1000;
cfg.method     = 'mlr';
% out3 = ft_denoise_tsr(cfg, data, refdata);
% 
% cfg.method     = 'mlrridge';
% cfg.threshold  = [0.1 0];
% out4 = ft_denoise_tsr(cfg, data, refdata);

cfg.testtrials = {1:10 11:20 21:30 31:40 41:50};
cfg.method     = 'mlr';
out5 = ft_denoise_tsr(cfg, data, refdata);

cfg.method     = 'mlrridge';
cfg.threshold  = [0.001 0];
out6 = ft_denoise_tsr(cfg, data, refdata);


x=1;
