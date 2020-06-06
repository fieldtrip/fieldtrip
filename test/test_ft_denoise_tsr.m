function test_ft_denoise_tsr

% MEM 8gb
% WALLTIME 01:30:00
% DEPENDENCY ft_denoise_tsr

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
cfg.reflags    = (-20:19)./1000;
cfg.method     = 'mlr';
out1 = ft_denoise_tsr(cfg, data, refdata);
assert(corr(out1.weights.beta(:),krn(:))>0.999 && norm(out1.weights.beta(:)-krn(:))./norm(krn(:))<0.01);

% now check the various demeaning/standardising stuff
cnt = 0;
for k1 = [true false]
  for k2 = [true false]
    for k3 = [true false]
      for k4 = [true false]
        cnt = cnt + 1;
        cfg.demeanrefdata      = k1;
        cfg.demeandata         = k2;
        cfg.standardiserefdata = k3;
        cfg.standardisedata    = k4;
        out1b(cnt) = ft_denoise_tsr(cfg, data, refdata);
      end
    end
  end
end
for k = 1:16
  w1b(k,:) = out1b(k).weights.beta;
end
assert(all(corr(w1b',krn(:))>0.999) && all(sqrt(sum((w1b-ones(16,1)*krn(:)').^2,2))<0.01));

cfg = [];
cfg.refchannel = {'refchan'};
cfg.reflags    = (-20:19)./1000;
cfg.method     = 'mlrridge';
cfg.standardiserefdata = true;
cfg.standardisedata    = true;
cfg.threshold  = [0.001 0];
out2 = ft_denoise_tsr(cfg, data, refdata);

data2 = data;
for k = 1:numel(data.trial)
  data2.trial{k} = data.trial{k} + randn(size(data.trial{k}));
end
cfg.threshold = [0.01 0];
out2n = ft_denoise_tsr(cfg, data2, refdata);

% cfg.method = 'mlrsvd';
% cfg.threshold = [1 1];
% out2svd = ft_denoise_tsr(cfg, data2, refdata);

for k = 1:numel(trial)
  trial{k} = trial{k} + 5.*ft_preproc_bandpassfilter(randn(1,1000), 1000, [8 12], [], 'firws');
end
data.trial = trial;

cfg = [];
cfg.refchannel = {'refchan'};
cfg.reflags    = (-20:19)./1000;
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
