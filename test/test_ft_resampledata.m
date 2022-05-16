function test_ft_resampledata

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_resampledata

fs = 500;
nchan = 32;
start_time = -1; % seconds
end_time = 2.5; % seconds
nsamples = (end_time - start_time) * fs + 1;

data = [];
data.time{1} = linspace(start_time, end_time, nsamples);
data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));

cfg = [];
cfg.resamplefs = 800;
dataout = ft_resampledata(cfg,data); % upsampling
cfg.resamplefs = 200;
dataout = ft_resampledata(cfg,data); % downsampling

% add a high amplitude broad-band component
data_hf_broad = data;
for k = 1:nchan
  data_hf_broad.trial{k}(1,:) = data_hf_broad.trial{1}(k,:) + ft_preproc_highpassfilter(randn(1,nsamples), fs, 100, [], 'firws').*10;
  data_hf_broad.time{k} = data_hf_broad.time{1};    
end
data_hf_broad.trial{1} = data_hf_broad.trial{1}(1,:);
data_hf_broad.label    = data_hf_broad.label(1);

% the default functionality in ft_resampledata applies a firls
% anti-aliasing filter that has its cutoff at the new Nyquist frequency
cfg = [];
cfg.resamplefs = 200;
dataout1 = ft_resampledata(cfg, data_hf_broad);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq   = 100;
dataout2 = ft_resampledata(cfg, data_hf_broad);

cfg.lpfreq = 90;
dataout3 = ft_resampledata(cfg, data_hf_broad);

cfg = [];
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.pad = 4;
freq = ft_freqanalysis(cfg, data_hf_broad);
freq1 = ft_freqanalysis(cfg, dataout1);
freq2 = ft_freqanalysis(cfg, dataout2);
freq3 = ft_freqanalysis(cfg, dataout3);

figure;hold on;
plot(freq.freq, (log10(freq.powspctrm)));
plot(freq1.freq, (log10(freq1.powspctrm)));
plot(freq2.freq, (log10(freq2.powspctrm)));
plot(freq3.freq, (log10(freq3.powspctrm)));

% add a high amplitude narrowband component
data_hf_narrow = data;
for k = 1:nchan
  data_hf_narrow.trial{k}(1,:) = data_hf_narrow.trial{1}(k,:) + ft_preproc_bandpassfilter(randn(1,nsamples), fs, [100 120], [], 'firws').*50;
  data_hf_narrow.time{k} = data_hf_narrow.time{1};    
end
data_hf_narrow.trial{1} = data_hf_narrow.trial{1}(1,:);
data_hf_narrow.label    = data_hf_narrow.label(1);

% the default functionality in ft_resampledata applies a firls
% anti-aliasing filter that has its cutoff at the new Nyquist frequency
cfg = [];
cfg.resamplefs = 200;
dataout1 = ft_resampledata(cfg, data_hf_narrow);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq   = 100;
dataout2 = ft_resampledata(cfg, data_hf_narrow);

cfg.lpfreq = 90;
dataout3 = ft_resampledata(cfg, data_hf_narrow);

cfg = [];
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.pad = 4;
freq = ft_freqanalysis(cfg, data_hf_narrow);
freq1 = ft_freqanalysis(cfg, dataout1);
freq2 = ft_freqanalysis(cfg, dataout2);
freq3 = ft_freqanalysis(cfg, dataout3);

figure;hold on;
plot(freq.freq, (log10(freq.powspctrm)));
plot(freq1.freq, (log10(freq1.powspctrm)));
plot(freq2.freq, (log10(freq2.powspctrm)));
plot(freq3.freq, (log10(freq3.powspctrm)));



