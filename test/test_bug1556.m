function test_bug1556

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_statfun_depsamplesFmultivariate ft_statfun_depsamplesFunivariate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1.powspctrm = [1     0     3     8     8    10     1     4     6     1]';
data1.dimord    = 'rpt_chan_freq_time';
data1.label     = {'cz'};
data1.freq      = 1;
data1.time      = 1;

data2.powspctrm = [1     8     6     1     6     8     7     4     9     0]';
data2.dimord    = 'rpt_chan_freq_time';
data2.label     = {'cz'};
data2.freq      = 1;
data2.time      = 1;

data3.powspctrm = [10     5     6     7    10     5    11     6    10     3]';
data3.dimord    = 'rpt_chan_freq_time';
data3.label     = {'cz'};
data3.freq      = 1;
data3.time      = 1;

cfg = [];
cfg.verbose  = 'off';
cfg.method   = 'analytic';
cfg.feedback = 'no';
cfg.alpha    = 5.0000e-02;
cfg.tail = 1;
cfg.statistic = 'ft_statfun_depsamplesFmultivariate';
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.design    = [ 1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2  3 3 3 3 3 3 3 3 3 3
  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 ];

stats = ft_freqstatistics(cfg, data1, data2, data3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.verbose  = 'off';
cfg.method   = 'analytic';
cfg.feedback = 'no';
cfg.alpha    = 5.0000e-02;
cfg.tail = 1;
cfg.design    = [ 1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2
                  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 ];

cfg.ivar      = 1;
cfg.statistic = 'ft_statfun_indepsamplesT';
stats_indepsamplesT = ft_freqstatistics(cfg, data1, data2);
cfg.statistic = 'ft_statfun_indepsamplesF';
stats_indepsamplesF = ft_freqstatistics(cfg, data1, data2);
cfg.uvar = 2;
cfg.statistic = 'ft_statfun_depsamplesT';
stats_depsamplesT = ft_freqstatistics(cfg, data1, data2);
cfg.statistic = 'ft_statfun_depsamplesFmultivariate';
stats_depsamplesF = ft_freqstatistics(cfg, data1, data2);

stats_indepsamplesT.stat
stats_indepsamplesF.stat
stats_depsamplesT.stat
stats_depsamplesF.stat

dat = cat(1, data1.powspctrm, data2.powspctrm);
x = data1.powspctrm;
y = data2.powspctrm;

[p, anovatab, stats_f] = anova1([x, y], []);
[h, p, ci, stats_t2] = ttest2(x, y);
[h, p, ci, stats_t] = ttest(x-y);

if abs(stats_t2.tstat^2/anovatab{2,5} - 1) > 0.001
  error('t^2 is unequal to F');
end


%% now testing the ft_statfun_depsamplesFunivariate

cfg = [];
cfg.verbose  = 'off';
cfg.method   = 'analytic';
cfg.feedback = 'no';
cfg.alpha    = 5.0000e-02;
cfg.tail = 1;
cfg.statistic = 'ft_statfun_depsamplesFunivariate';% computations are based on sums of squares formula
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.design    = [ 1 1 1 1 1 1 1 1 1 1  2 2 2 2 2 2 2 2 2 2  3 3 3 3 3 3 3 3 3 3
                  1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 ];

stats = ft_freqstatistics(cfg, data1, data2, data3);



%% Contrasting ANOVA by GLM. Should be equivalent to the sums of squares way 
% Source: visit http://www.sbirc.ed.ac.uk/cyril/glm/GLM_lectures.html#10
% for an excellent tutorial on GLM and here to find formulas
% http://www.fil.ion.ucl.ac.uk/~wpenny/publications/rik_anova.pdf

Y = [data1.powspctrm; data2.powspctrm; data3.powspctrm];

% create the design matrix for the different factors
nb_subjects =10;
nb_conditions =3;
Subjects = repmat(eye(nb_subjects),nb_conditions,1); % error
x = kron(eye(nb_conditions),ones(nb_subjects,1));  % effect
X = [x Subjects]; % no more ones for the grand mean but a subject specific mean
figure; imagesc(X); colormap('gray'); title('Repearted measure design','Fontsize',14)

% Compute as usual
df  = nb_conditions -1;
dfe = size(Y,1)  - nb_subjects - df;

P     = X*pinv(X'*X)*X'; % our projection matrix
R     = eye(size(Y,1)) - P; % projection on error space
SSe   = diag(Y'*R*Y); % Y projected onto the error
Betas = pinv(x)*Y;  % compute without cst/subjects
yhat  = x*Betas; % yhat computed based on the treatment with subject effect - we use little x
SS    = norm(yhat-mean(yhat)).^2;
F_values = (SS/df) ./ (SSe/dfe);
p_values = 1 - fcdf(F_values, df, dfe);


%% now comparing the two appraches
tolerance = 1e-4;
if abs(stats.stat - F_values) > tolerance;
  error('F univariate statistics are unequal');
end
if abs(stats.prob - p_values) > tolerance;
  error('pvalues differ');
end
if abs(stats.dfdenom - dfe) > tolerance;
  error('degrees of freedom of the error term differ');
end
if abs(stats.dfnum - df) > tolerance;
  error('degrees of freedom of the factor term differ');
end
