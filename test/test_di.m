%
% WALLTIME:
% MEM: 
% 
% 

%% Simulate the signals (10 channels)

fs = 150;           % in Hz
total = 3;          % in seconds
t = (0:1/fs:total); % time axis

% connectivity timecourse
starttime = 0.2;  % when source arrives
offtime = 0.3;    % when source is off

% preallocate data matrix
sim = zeros(10, size(t, 2));

% define signal components

omega = 0.1;
drift = 0.1*sin(2*pi*t*omega); 
omega = 4;
theta = 0.4*sin(2*pi*t*omega);
omega = 40;
gamma = 0.1*sin(2*pi*t*omega);

% Create 2 source signals with noisepl
seeds = [123, 321];
sources = zeros(2, length(t));
for i = 1:numel(seeds)
    
    rng(seeds(i));
    %noise1 = smoothdata(randn(1, length(t)).*(0.2.*-1.*hamming(length(t))+1)', 'gaussian', 30);
    noise1 = ft_preproc_smooth(randn(1, length(t)).*(0.2.*-1.*hamming(length(t))+1)', 30);
    
    noise2 = 0.2*randn(1, length(t));

    sources(i, :) = drift + theta.*hamming(length(t))' + gamma + noise1 + noise2;

end

% fill unrelated noise signals
rng(100);
tap = repmat((0.2.*-1.*hamming(length(t))+1)', [5, 1]);
%noise1tmp = smoothdata((randn(5, length(t)).*tap), 2, 'gaussian', 30);
noise1tmp = ft_preproc_smooth((randn(5, length(t)).*tap), 30);
noise2tmp = (0.2*randn(5, length(t)));
sim(1:2, :) = noise1tmp(1:2, :) + noise2tmp(1:2, :);
sim(5:7, :) = noise1tmp(3:5, :) + noise2tmp(3:5, :);

% fill sources 
sim(3, :) = sources(1, :);
sim(4, :) = sources(2, :);

% define the start and end sample points
endsmp = nearest(t, 10);
begsmp = nearest(t, 4);

% create targets by making them equal to the past of the sources
targets = zeros(3, endsmp);
targets(1, :) = sources(1, 1:endsmp);
targets(2, :) = sources(2, 1:endsmp);
targets(3, :) = 0.5.*sources(1, 1:endsmp) + 0.5.*sources(2, 1:endsmp);

sim = sim(:, 1:endsmp);

% write targets to the data matrix
sim(8, :) = targets(1, :);
sim(9, :) = targets(2, :);
sim(10, :) = targets(3, :);

%% Plot the signals

tnew = t(1:endsmp);

figure();
subplot(4, 1, 1);
plot(tnew, sim(3, :));
legend('source 1');
subplot(4, 1, 2);
plot(tnew, sim(4, :));
subplot(4, 1, 3);
legend('source 2');
plot(tnew, sim(8, :));
line([tnew(begsmp), tnew(begsmp)], [-1, 1], 'color', 'red');
legend('target 1', 'source 1 arrives');
subplot(4, 1, 4);
plot(tnew, sim(9, :));
line([tnew(begsmp), tnew(begsmp)], [-1, 1], 'color', 'red');
legend('target 2', 'source 2 arrives');
xlabel('time (sec)');


% plot the whole matrix ('sim' variable)

figure();
imagesc(sim);
ax = gca;
ax.YTickLabel{3} = 'S1';
ax.YTickLabel{4} = 'S2';
ax.YTickLabel{8} = 'S1-tgt';
ax.YTickLabel{9} = 'S2-tgt';
ax.YTickLabel{10} = 'S1+S2-tgt';
title('Data');
ylabel('channels');
xlabel('time (sample point)');


%% create labels

labels = cell(size(sim, 1),1);
for j = 1:size(sim, 1)
    labels{j} = sprintf('chan%d', j); 
end

% create ft-style struct

data.fsample = fs;
data.trial = {sim+randn(size(sim))./100};
data.time = {t(1:endsmp)};
data.label = labels;
data.sampleinfo = [1, endsmp];

%% Run 
cfg            = [];
cfg.method     = 'di';
cfg.refindx    = 'all';
cfg.di.lags    = (0.05:0.05:0.3);

% di = ft_connectivityanalysis(cfg, data);
% 
% %% Plot connectivity time courses
% 
% figure;
% title('Di')
% plot(di.time, squeeze(di.di(3, 8, :)), '-o'); hold on;
% plot(di.time, squeeze(di.di(9, 4, :)), '-o');
% plot(di.time, squeeze(di.di(8, 3, :)), '-o');
% plot(di.time, squeeze(di.di(3, 6, :)), '-o');
% ylabel('di (bit)')
% xlabel('lag (sec)');
% labs = {'3(source)-->8(target)', '4(source)-->9(target)', '8(target)-->3(source)', '3(source)-->6(non-target)'};
% legend(labs);
