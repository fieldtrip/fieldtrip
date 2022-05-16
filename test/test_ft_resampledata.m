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

% convert the data into a single channel, multi trial representation so
% that a meaningful powerspectrum can be computed later on
data_ = data;
for k = 1:nchan
  data_.trial{k}(1,:) = data_.trial{1}(k,:);
  data_.time{k} = data_.time{1};    
end
data_.trial{1} = data_.trial{1}(1,:);
data_.label    = data_.label(1);

% the default functionality in ft_resampledata applies a firls
% anti-aliasing filter that has its cutoff at the new Nyquist frequency
cfg = [];
cfg.resamplefs = 200;
dataout1 = ft_resampledata(cfg, data_);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq   = 100;
dataout2 = ft_resampledata(cfg, data_);

cfg.lpfreq = 90;
dataout3 = ft_resampledata(cfg, data_);

cfg = [];
cfg.time = dataout1.time;
dataout4 = ft_resampledata(cfg, data_);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout5 = ft_resampledata(cfg, data_);

cfg = [];
cfg.resamplefs = 250;
cfg.method = 'downsample';
dataout6 = ft_resampledata(cfg, data_);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout7 = ft_resampledata(cfg, data_);

cfg = [];
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.pad = 4;
freq = ft_freqanalysis(cfg, data_);
freq1 = ft_freqanalysis(cfg, dataout1);
freq2 = ft_freqanalysis(cfg, dataout2);
freq3 = ft_freqanalysis(cfg, dataout3);
freq4 = ft_freqanalysis(cfg, dataout4);
freq5 = ft_freqanalysis(cfg, dataout5);
freq6 = ft_freqanalysis(cfg, dataout6);
freq7 = ft_freqanalysis(cfg, dataout7);

cmap = ft_colormap('Set1');

figure;hold on;
plot(freq1.freq, (log10(freq1.powspctrm)), 'color', cmap(2,:), 'linewidth', 2);
plot(freq2.freq, (log10(freq2.powspctrm)), 'color', cmap(3,:), 'linewidth', 2);
plot(freq3.freq, (log10(freq3.powspctrm)), 'color', cmap(4,:), 'linewidth', 2);
plot(freq4.freq, (log10(freq4.powspctrm)), 'color', cmap(5,:), 'linewidth', 2);
plot(freq5.freq, (log10(freq5.powspctrm)), 'color', cmap(7,:), 'linewidth', 2);
plot(freq6.freq, (log10(freq6.powspctrm)), 'color', cmap(8,:), 'linewidth', 2);
plot(freq7.freq, (log10(freq7.powspctrm)), 'color', cmap(9,:), 'linewidth', 2);
plot(freq.freq,  (log10(freq.powspctrm)), 'color', cmap(1,:), 'linewidth', 2);
legend({'rs_native', 'rs_firws100', 'rs_firws090', 'interp1', 'interp1_firws100', 'downsample', 'downsample_firws100', 'original'}, 'interpreter', 'none');
xlabel('frequency (Hz)');
ylabel('power');

% add a high amplitude broad-band component
data_hf_broad = data;
for k = 1:nchan
  data_hf_broad.trial{k}(1,:) = data_hf_broad.trial{1}(k,:) + ft_preproc_highpassfilter(randn(1,nsamples), fs, 110, [], 'firws').*50;
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
cfg.time = dataout1.time;
dataout4 = ft_resampledata(cfg, data_hf_broad);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout5 = ft_resampledata(cfg, data_hf_broad);

cfg = [];
cfg.resamplefs = 250;
cfg.method = 'downsample';
dataout6 = ft_resampledata(cfg, data_hf_broad);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout7 = ft_resampledata(cfg, data_hf_broad);

cfg = [];
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.pad = 4;
freq = ft_freqanalysis(cfg, data_hf_broad);
freq1 = ft_freqanalysis(cfg, dataout1);
freq2 = ft_freqanalysis(cfg, dataout2);
freq3 = ft_freqanalysis(cfg, dataout3);
freq4 = ft_freqanalysis(cfg, dataout4);
freq5 = ft_freqanalysis(cfg, dataout5);
freq6 = ft_freqanalysis(cfg, dataout6);
freq7 = ft_freqanalysis(cfg, dataout7);

figure;hold on;
plot(freq1.freq, (log10(freq1.powspctrm)), 'color', cmap(2,:), 'linewidth', 2);
plot(freq2.freq, (log10(freq2.powspctrm)), 'color', cmap(3,:), 'linewidth', 2);
plot(freq3.freq, (log10(freq3.powspctrm)), 'color', cmap(4,:), 'linewidth', 2);
plot(freq4.freq, (log10(freq4.powspctrm)), 'color', cmap(5,:), 'linewidth', 2);
plot(freq5.freq, (log10(freq5.powspctrm)), 'color', cmap(7,:), 'linewidth', 2);
plot(freq6.freq, (log10(freq6.powspctrm)), 'color', cmap(8,:), 'linewidth', 2);
plot(freq7.freq, (log10(freq7.powspctrm)), 'color', cmap(9,:), 'linewidth', 2);
plot(freq.freq,  (log10(freq.powspctrm)), 'color', cmap(1,:), 'linewidth', 2);
legend({'rs_native', 'rs_firws100', 'rs_firws090', 'interp1', 'interp1_firws100', 'downsample', 'downsample_firws100', 'original'}, 'interpreter', 'none');
xlabel('frequency (Hz)');
ylabel('power');

% add a high amplitude narrowband component
data_hf_narrow = data;
for k = 1:nchan
  data_hf_narrow.trial{k}(1,:) = data_hf_narrow.trial{1}(k,:) + ft_preproc_bandpassfilter(randn(1,nsamples), fs, [110 120], [], 'firws').*50;
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
cfg.time = dataout1.time;
dataout4 = ft_resampledata(cfg, data_hf_narrow);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout5 = ft_resampledata(cfg, data_hf_narrow);

cfg = [];
cfg.resamplefs = 250;
cfg.method = 'downsample';
dataout6 = ft_resampledata(cfg, data_hf_narrow);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout7 = ft_resampledata(cfg, data_hf_broad);

cfg = [];
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.pad = 4;
freq = ft_freqanalysis(cfg, data_hf_narrow);
freq1 = ft_freqanalysis(cfg, dataout1);
freq2 = ft_freqanalysis(cfg, dataout2);
freq3 = ft_freqanalysis(cfg, dataout3);
freq4 = ft_freqanalysis(cfg, dataout4);
freq5 = ft_freqanalysis(cfg, dataout5);
freq6 = ft_freqanalysis(cfg, dataout6);
freq7 = ft_freqanalysis(cfg, dataout7);

figure;hold on;
plot(freq1.freq, (log10(freq1.powspctrm)), 'color', cmap(2,:), 'linewidth', 2);
plot(freq2.freq, (log10(freq2.powspctrm)), 'color', cmap(3,:), 'linewidth', 2);
plot(freq3.freq, (log10(freq3.powspctrm)), 'color', cmap(4,:), 'linewidth', 2);
plot(freq4.freq, (log10(freq4.powspctrm)), 'color', cmap(5,:), 'linewidth', 2);
plot(freq5.freq, (log10(freq5.powspctrm)), 'color', cmap(7,:), 'linewidth', 2);
plot(freq6.freq, (log10(freq6.powspctrm)), 'color', cmap(8,:), 'linewidth', 2);
plot(freq7.freq, (log10(freq7.powspctrm)), 'color', cmap(9,:), 'linewidth', 2);
plot(freq.freq,  (log10(freq.powspctrm)), 'color', cmap(1,:), 'linewidth', 2);
legend({'rs_native', 'rs_firws100', 'rs_firws090', 'interp1', 'interp1_firws100', 'downsample', 'downsample_firws100', 'original'}, 'interpreter', 'none');
xlabel('frequency (Hz)');
ylabel('power');


% make 1/f
data_1f = data;
for k = 1:nchan
  data_1f.trial{k}(1,:) = data_1f.trial{1}(k,:) + cumsum(randn(1,nsamples))./10;
  data_1f.trial{k} = data_1f.trial{k} - mean(data_1f.trial{k},2);
  data_1f.time{k} = data_1f.time{1};    
end
data_1f.trial{1} = data_1f.trial{1}(1,:);
data_1f.label    = data_1f.label(1);

% the default functionality in ft_resampledata applies a firls
% anti-aliasing filter that has its cutoff at the new Nyquist frequency
cfg = [];
cfg.resamplefs = 200;
dataout1 = ft_resampledata(cfg, data_1f);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq   = 100;
dataout2 = ft_resampledata(cfg, data_1f);

cfg.lpfreq = 90;
dataout3 = ft_resampledata(cfg, data_1f);

cfg = [];
cfg.time = dataout1.time;
dataout4 = ft_resampledata(cfg, data_1f);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout5 = ft_resampledata(cfg, data_1f);

cfg = [];
cfg.resamplefs = 250;
cfg.method = 'downsample';
dataout6 = ft_resampledata(cfg, data_1f);

cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 100;
dataout7 = ft_resampledata(cfg, data_1f);

cfg = [];
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.pad = 4;
freq = ft_freqanalysis(cfg, data_1f);
freq1 = ft_freqanalysis(cfg, dataout1);
freq2 = ft_freqanalysis(cfg, dataout2);
freq3 = ft_freqanalysis(cfg, dataout3);
freq4 = ft_freqanalysis(cfg, dataout4);
freq5 = ft_freqanalysis(cfg, dataout5);
freq6 = ft_freqanalysis(cfg, dataout6);
freq7 = ft_freqanalysis(cfg, dataout7);

figure;hold on;
plot(freq1.freq, (log10(freq1.powspctrm)), 'color', cmap(2,:), 'linewidth', 2);
plot(freq2.freq, (log10(freq2.powspctrm)), 'color', cmap(3,:), 'linewidth', 2);
plot(freq3.freq, (log10(freq3.powspctrm)), 'color', cmap(4,:), 'linewidth', 2);
plot(freq4.freq, (log10(freq4.powspctrm)), 'color', cmap(5,:), 'linewidth', 2);
plot(freq5.freq, (log10(freq5.powspctrm)), 'color', cmap(7,:), 'linewidth', 2);
plot(freq6.freq, (log10(freq6.powspctrm)), 'color', cmap(8,:), 'linewidth', 2);
plot(freq7.freq, (log10(freq7.powspctrm)), 'color', cmap(9,:), 'linewidth', 2);
plot(freq.freq,  (log10(freq.powspctrm)), 'color', cmap(1,:), 'linewidth', 2);
legend({'rs_native', 'rs_firws100', 'rs_firws090', 'interp1', 'interp1_firws100', 'downsample', 'downsample_firws100', 'original'}, 'interpreter', 'none');
xlabel('frequency (Hz)');
ylabel('power');
