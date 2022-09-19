function test_bug3417

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY

% ... works fine on his computer, which has a 2012 version of FieldTrip and Matlab 2013a.
% But when we try to run either a later version of Fieldtrip (e.g. 2015) or a later version
% of Matlab (e.g. 2017b), we get the following bug:
%
% Error using findcluster (line 50)
% invalid dimension of spatdimneighbstructmat

% figure 4

%% SL - one sample t

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3417.mat'))

freqs=14; % or 21
times=81;

cond1=permute(SL_Fig4,[3,1,2]);
cond1=cond1(:,1:freqs,:);

cond2=cond1;
cond2(:,:,:)=0;

data1.globspctrm=cond1;
data2.globspctrm=cond2;

[subj,~,freq] = size(data1.globspctrm);
nsubj = subj;

cfg.statistic = 'depsamplesT';
cfg.design = zeros(2,2*nsubj);
cfg.design(1,:)=[ones(1,nsubj) ones(1,nsubj)+1];
cfg.design(2,:)=[1:nsubj 1:nsubj];
cfg.ivar = 1;
cfg.uvar = 2;

data1.label = {'Positive'};
data2.label = {'Negative'};
data1.dimord = 'rpt_freq_time'; % does clustering in 2 dimensions (freq, time)
data2.dimord = 'rpt_freq_time';
data1.freq = [1:freqs];
data2.freq = [1:freqs];
data1.time =[1:times];
data2.time= [1:times];

% fix the invalid data structure, see ft_datatype_freq
data1.dimord = 'rpt_chan_freq_time';
data2.dimord = 'rpt_chan_freq_time';
data1.globspctrm = reshape(data1.globspctrm, [nsubj 1 freqs times]);
data2.globspctrm = reshape(data2.globspctrm, [nsubj 1 freqs times]);
data1.label = {'TheSame'};
data2.label = {'TheSame'};

cfg.parameter = 'globspctrm'; % 'avg', 'cohspctrm'
cfg.method = 'montecarlo'; %  'montecarlo',   'analytic' , 'stats', 'crossvalidate'

cfg.neighbours = [];
cfg.channel = 'all';
cfg.correctm = 'cluster'; % 'no', 'max', 'cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')

cfg.clusterstatistc = 'wcm'; % 'wcm', 'maxsum', 'maxsize'
cfg.clusterthreshold = 'parametric'; %  'parametric', 'nonparametric_individual', 'nonparametric_common'
cfg.numrandomization = 1000;
cfg.clusteralpha = 0.05;
cfg.tail = 0;

stats = ft_freqstatistics(cfg,data1,data2);

figure;

subplot(2,3,1)
tres=squeeze(mean(cond1,1));
imagesc(tres)
colorbar
title('Adult SL: Mean')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,4)
tres=squeeze(stats.stat)
imagesc(tres)
colorbar
title('SL: t matrix')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,5)
Positive_cluster=squeeze(stats.posclusterslabelmat);
Positive_cluster(Positive_cluster > 1) = 0; % keep just cluster 1
imagesc(Positive_cluster)
title(['Positive cluster: p=' num2str(stats.posclusters(1, 1).prob)])
colorbar
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,6)
Negative_cluster=squeeze(stats.negclusterslabelmat);
Negative_cluster(Negative_cluster > 1) = 0; % keep just cluster 1
imagesc(Negative_cluster)
title(['Negative cluster: p=' num2str(stats.negclusters(1, 1).prob)])
colorbar
ylabel('Frequencies')
xlabel('Time')


%% JA - one sample t

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3417.mat'))

freqs=14; % or 21
times=81;

cond1=permute(JA_Fig4,[3,1,2]);
cond1=cond1(:,1:freqs,:);

data1.globspctrm=cond1;
data2.globspctrm=cond2;


[subj,~,freq] = size(data1.globspctrm);
nsubj = subj;

cfg.statistic = 'depsamplesT';
cfg.design = zeros(2,2*nsubj);
cfg.design(1,:)=[ones(1,nsubj) ones(1,nsubj)+1];
cfg.design(2,:)=[1:nsubj 1:nsubj];
cfg.ivar = 1;
cfg.uvar = 2;

data1.label = {'Positive'};
data2.label = {'Negative'};
data1.dimord = 'rpt_freq_time'; % does clustering in 2 dimensions (freq, time)
data2.dimord = 'rpt_freq_time';
data1.freq = [1:freqs];
data2.freq = [1:freqs];
data1.time =[1:times];
data2.time= [1:times];

% fix the invalid data structure, see ft_datatype_freq
data1.dimord = 'rpt_chan_freq_time';
data2.dimord = 'rpt_chan_freq_time';
data1.globspctrm = reshape(data1.globspctrm, [nsubj 1 freqs times]);
data2.globspctrm = reshape(data2.globspctrm, [nsubj 1 freqs times]);
data1.label = {'TheSame'};
data2.label = {'TheSame'};

cfg.parameter = 'globspctrm'; % 'avg', 'cohspctrm'
cfg.method = 'montecarlo'; %  'montecarlo',   'analytic' , 'stats', 'crossvalidate'

cfg.neighbours = [];
cfg.channel = 'all';
cfg.correctm = 'cluster'; % 'no', 'max', 'cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')

cfg.clusterstatistc = 'wcm'; % 'wcm', 'maxsum', 'maxsize'
cfg.clusterthreshold = 'parametric'; %  'parametric', 'nonparametric_individual', 'nonparametric_common'
cfg.numrandomization = 1000;
cfg.clusteralpha = 0.05;
cfg.tail = 0;

stats = ft_freqstatistics(cfg,data1,data2)

figure;

subplot(2,3,1)
tres=squeeze(mean(cond1,1));
imagesc(tres)
colorbar
title('Adult SL: Mean')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,4)
tres=squeeze(stats.stat)
imagesc(tres)
colorbar
title('SL: t matrix')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,5)
Positive_cluster=squeeze(stats.posclusterslabelmat);
Positive_cluster(Positive_cluster > 1) = 0; % keep just cluster 1
imagesc(Positive_cluster)
title(['Positive cluster: p=' num2str(stats.posclusters(1, 1).prob)])
colorbar
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,6)
Negative_cluster=squeeze(stats.negclusterslabelmat);
Negative_cluster(Negative_cluster > 1) = 0; % keep just cluster 1
imagesc(Negative_cluster)
title(['Negative cluster: p=' num2str(stats.negclusters(1, 1).prob)])
colorbar
ylabel('Frequencies')
xlabel('Time')

%% SJ vs. AJ

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3417.mat'))

freqs=14; % or 21
times=81;

cond1=permute(SL_Fig4,[3,1,2]);
cond2=permute(JA_Fig4,[3,1,2]);

cond1=cond1(:,1:freqs,:);
cond2=cond2(:,1:freqs,:);

data1.globspctrm=cond1;
data2.globspctrm=cond2;


[subj,~,freq] = size(data1.globspctrm);
nsubj = subj;

cfg.statistic = 'depsamplesT';
cfg.design = zeros(2,2*nsubj);
cfg.design(1,:)=[ones(1,nsubj) ones(1,nsubj)+1];
cfg.design(2,:)=[1:nsubj 1:nsubj];
cfg.ivar = 1;
cfg.uvar = 2;

data1.label = {'Positive'};
data2.label = {'Negative'};
data1.dimord = 'rpt_freq_time'; % does clustering in 2 dimensions (freq, time)
data2.dimord = 'rpt_freq_time';
data1.freq = [1:freqs];
data2.freq = [1:freqs];
data1.time =[1:times];
data2.time= [1:times];

% fix the invalid data structure, see ft_datatype_freq
data1.dimord = 'rpt_chan_freq_time';
data2.dimord = 'rpt_chan_freq_time';
data1.globspctrm = reshape(data1.globspctrm, [nsubj 1 freqs times]);
data2.globspctrm = reshape(data2.globspctrm, [nsubj 1 freqs times]);
data1.label = {'TheSame'};
data2.label = {'TheSame'};

cfg.parameter = 'globspctrm'; % 'avg', 'cohspctrm'
cfg.method = 'montecarlo'; %  'montecarlo',   'analytic' , 'stats', 'crossvalidate'

cfg.neighbours = [];
cfg.channel = 'all';
cfg.correctm = 'cluster'; % 'no', 'max', 'cluster', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')

cfg.clusterstatistc = 'wcm'; % 'wcm', 'maxsum', 'maxsize'
cfg.clusterthreshold = 'parametric'; %  'parametric', 'nonparametric_individual', 'nonparametric_common'
cfg.numrandomization = 1000;
cfg.clusteralpha = 0.05;
cfg.tail = 0;

stats = ft_freqstatistics(cfg,data1,data2)

figure;

subplot(2,3,1)
tres=squeeze(mean(cond1,1));
imagesc(tres)
colorbar
title('Adult SL: Mean')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,2)
tres=squeeze(mean(cond2,1));
imagesc(tres)
colorbar
title('Adult JA: Mean')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,4)
tres=squeeze(stats.stat)
imagesc(tres)
colorbar
title('SL vs. JA: t matrix')
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,5)
Positive_cluster=squeeze(stats.posclusterslabelmat);
Positive_cluster(Positive_cluster > 1) = 0; % keep just cluster 1
imagesc(Positive_cluster)
title(['SL > JA cluster: p=' num2str(stats.posclusters(1, 1).prob)])
colorbar
ylabel('Frequencies')
xlabel('Time')

subplot(2,3,6)
Negative_cluster=squeeze(stats.negclusterslabelmat);
Negative_cluster(Negative_cluster > 1) = 0; % keep just cluster 1
imagesc(Negative_cluster)
title(['SL < JA cluster: p=' num2str(stats.negclusters(1, 1).prob)])
colorbar
ylabel('Frequencies')
xlabel('Time')
