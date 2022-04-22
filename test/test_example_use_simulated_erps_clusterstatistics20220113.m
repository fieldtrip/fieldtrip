function test_example_use_simulated_erps_to_explore_cluster_statistics

% MEM 4gb
% WALLTIME 00:10:00

%
%% Use simulated ERPs to explore cluster statistics
%
% The following code starts off with an ERP in two conditions, where it is slightly larger in condition 1 than 2. This simulation demonstrates a randomization test, correcting for multiple comparisons by using the largest cluster mass.
%
% See this paper for more details
%
% Maris E., Oostenveld R. [Nonparametric statistical testing of EEG- and MEG-data.](http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=17517438) J Neurosci Methods. 2007 Apr 10;
%
% and look in the [reference](/references_to_implemented_methods) section for more literature pointers.
%
%%

tim = (0:1000)/1000;
erp = (1-cos(2*pi*tim))/2;
plot(tim, erp)
xlabel('time (s)')
ylabel('ERP (uV)')

%%

data1 = [];
for i=1:100
  data1.time{i} = tim;
  data1.trial{i}(1,:) = 1.1*erp + 0.1*randn(size(erp)); % the effect size is specified here
end
data1.label = {'Cz'};

%
data2 = [];
for i=1:100
  data2.time{i} = tim;
  data2.trial{i}(1,:) = 1.0*erp + 0.1*randn(size(erp));
end
data2.label = {'Cz'};

%%

cfg = [];
cfg.keeptrials = 'no';
avg1 = ft_timelockanalysis(cfg, data1);
avg2 = ft_timelockanalysis(cfg, data2);

figure
cfg = [];
ft_singleplotER(cfg, avg1, avg2);
legend({'condition 1', 'condition 2'})

%%

cfg = [];
cfg.operation = 'x1-x2';
cfg.parameter = 'avg';
difference = ft_math(cfg, avg1, avg2);

figure
cfg = [];
ft_singleplotER(cfg, difference);
legend({'difference'})

%%

% for the statistics we need the variance over the individual trials
cfg = [];
cfg.keeptrials = 'yes';
timelock1 = ft_timelockanalysis(cfg, data1);
timelock2 = ft_timelockanalysis(cfg, data2);

%
cfg = [];
cfg.design = [ 1*ones(1,100) 2*ones(1,100) ];
cfg.ivar = 1;
cfg.method = 'montecarlo';
cfg.statistic = 'indepsamplesT';
cfg.correctm = 'cluster';
cfg.numrandomization = 2000;
% cfg.neighbours = []; % only cluster over time, not over channels
cfg.spmversion = 'spm12'; % the default spm8 has mex file problems on recent macOS versions
stat = ft_timelockstatistics(cfg, timelock1, timelock2);

%%

figure
subplot(4,1,1); plot(stat.time, stat.stat); ylabel('t-value');
subplot(4,1,2); plot(difference.time, difference.avg); ylabel('avg1-avg2 (uV)');
subplot(4,1,3); semilogy(stat.time, stat.prob); ylabel('prob'); axis([0 1 0.001 2])
subplot(4,1,4); plot(stat.time, stat.mask); ylabel('significant'); axis([0 1 -0.1 1.1])

%%

figure
subplot(2,1,1); hist(stat.negdistribution, 200); axis([-10 10 0 100])
for i=1:numel(stat.negclusters)
  X = [stat.negclusters(i).clusterstat stat.negclusters(i).clusterstat];
  Y = [0 100];
  line(X, Y, 'color', 'r')
end

subplot(2,1,2); hist(stat.posdistribution, 200); axis([-10 10 0 100])

for i=1:numel(stat.posclusters)
  X = [stat.posclusters(i).clusterstat stat.posclusters(i).clusterstat];
  Y = [0 100];
  line(X, Y, 'color', 'r')
end
