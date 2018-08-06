function test_bug2364

% WALLTIME 00:10:00
% MEM 1gb


nTrials = 20;
data = [];
data.fsample = 256;
data.trial = arrayfun(@(x) rand(10, 512), 1:nTrials, 'UniformOutput', false);
data.time = repmat({(0:length(data.trial{1})-1)/data.fsample}, [1, nTrials]);
data.label = cellfun(@num2str, num2cell(1:10), 'UniformOutput', false);

for iTrial = 1:nTrials
  data.trial{iTrial}(1,:) = nan;
end


%%

cfg = [];
cfg.method = 'wavelet';
cfg.width = 2;
cfg.foilim = [10 50];
cfg.toi = 0:0.01:2;

freq1 = ft_freqanalysis(cfg, data);
count1 = sum(~isnan(freq1.powspctrm(:)));
fprintf('number of non-nans for wavelet method: %i\n', count1)


%%

cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 10:50;
cfg.t_ftimwin = 2./cfg.foi;
cfg.toi = 0:0.01:2;

freq2 = ft_freqanalysis(cfg, data);
count2 = sum(~isnan(freq2.powspctrm(:)));
fprintf('number of non-nans for mtmconvol method: %i\n', count2)

%%

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';

freq3 = ft_freqanalysis(cfg, data);
count3 = sum(~isnan(freq3.powspctrm(:)));
fprintf('number of non-nans for mtmfft method: %i\n', count3)

%%

assert(count1==140067)
assert(count2==70803)
assert(count3==2313);


%% deal with regression problem in single-channel data

data1 = [];
data1.label = {'1'};
data1.time = {(1:1000)/1000};
data1.trial = {randn(1,1000)};

data2 = [];
data2.label = {'1', '2'};
data2.time = {(1:1000)/1000};
data2.trial = {randn(2,1000)};
data2.trial{1}(1,:) = data1.trial{1};

cfg = [];
cfg.method = 'wavelet';
cfg.width = 2;
cfg.foilim = [10 50];
cfg.toi = 0:0.01:1;
tf1 = ft_freqanalysis(cfg, data1);
tf2 = ft_freqanalysis(cfg, data2);

figure; imagesc(abs(squeeze(tf1.powspctrm(1,:,:))))
figure; imagesc(abs(squeeze(tf2.powspctrm(1,:,:))))

difference = abs(tf1.powspctrm(1,:,:) - tf2.powspctrm(1,:,:));
assert(all(difference(~isnan(difference(:)))<100*eps)); % there is a small numerical difference



