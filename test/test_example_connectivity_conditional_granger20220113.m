function test_example_conditional_granger

% MEM 4gb
% WALLTIME 00:10:00

%
%% Conditional Granger causality in the frequency domain
%
% Conditional Granger causality is a derivative of spectral Granger causality that is computed over a triplet of channels (or blocks of channels). It provides the advantage that for this triplet, it allows to differentiate between a delayed parallel drive from sources <i>A</i> to be <i>B</i> and <i>C</i> and a sequential drive from <i>A</i> to <i>B</i> to <i>C</i>.
%
% This example illustrates the simulation and base analysis of the paper
%
% See also: [the connectivity tutorial]/tutorial/connectivity/).
%
%% # Setup and simulating the data sets
%
% First, define parameters under which samples should be simulated.
%
simcfg             = [];
simcfg.ntrials     = 500;
simcfg.triallength = 1;
simcfg.fsample     = 200;
simcfg.nsignal     = 3;
simcfg.method      = 'ar';

% We want to simulate a system with three signals. Their noise is modeled as white noise processes with zero mean and standard deviations
%
%
% They will have the covariances &zeta;, &eta; and &epsilon;. We also require a paramters &mu;=0.5.
%
% parameters of the model itself
mu                 = 0.5;
absnoise           = [ 1.0   0.2   0.3 ];

% First, we generate the sample for the case of sequential driving.
% We want to incorporate the system
%
%
% which we can do like this:
%
% params(i,j,k): j -> i at t=k
simcfg.params(:,:,1) = [   0      0      0;
                         1.0      0      0;
                           0    1.0     mu];

% Note that the matrix representation for the covariance reads from columns to row, other than the MVAR-model is read intuitively. But we still need to hand the parameters of the noise to the model:
%
% paper defines stds, not cov:
simcfg.noisecov      = diag(absnoise.^2);

data2           = ft_connectivitysimulation(simcfg);

% Now create sample data for the case of differentially delayed driving,
%
%
% which we can write as
%
simcfg.params(:,:,1) = [   0      0      0;
                         1.0      0      0;
                           0      0     mu];
simcfg.params(:,:,2) = [   0      0      0;
                           0      0      0;
                         1.0      0      0];

% We build the actual MVAR-representation...
%
data1           = ft_connectivitysimulation(simcfg);

% #
%
figure
plot(data1.time{1}, data1.trial{1})
legend(data1.label)
xlabel('time (s)')

% Don't be confused that we started with data2 and conclude with data1. This is just to maintain the order the systems have in the paper.
%
%% # MVAR model frequency analysis
% We generate spectral representations from the MVAR representations we defined with data1 and data2. After all, we want to compute spectral Granger causality. Fast Fourier is a good starting point.
%
freq                   = [];
freq.freqcfg           = [];
freq.freqcfg.method    = 'mtmfft';
freq.freqcfg.output    = 'fourier';
freq.freqcfg.tapsmofrq = 2;
freqdata1           = ft_freqanalysis(freq.freqcfg, data1);
freqdata2           = ft_freqanalysis(freq.freqcfg, data2);

%% # "Regular" Granger causality
% Let first compute regular bivariate Granger causality, as this makes the difference clear to what we want.
%
grangercfg = [];
grangercfg.method  = 'granger';
grangercfg.granger.conditional = 'no';
grangercfg.granger.sfmethod = 'bivariate';

gdata = [];
gdata.g1_bivar_reg      = ft_connectivityanalysis(grangercfg, freqdata1);
gdata.g2_bivar_reg      = ft_connectivityanalysis(grangercfg, freqdata2);

%% # Multivariate conditional Granger causality
% However, we clearly want a multivariate approach. Also, we need to define channel combinations, as we now require triplets of inputs.
%
grangercfg.granger.conditional = 'yes';
grangercfg.channelcmb  = {'signal001', 'signal002', 'signal003'};
grangercfg.granger.sfmethod = 'multivariate';
grangercfg.granger.conditional = 'yes';

% block-wise causality
grangercfg.granger.block(1).name   = freqdata1.label{1};
grangercfg.granger.block(1).label  = freqdata1.label(1);
grangercfg.granger.block(2).name   = freqdata1.label{2};
grangercfg.granger.block(2).label  = freqdata1.label(2);
grangercfg.granger.block(3).name   = freqdata1.label{3};
grangercfg.granger.block(3).label  = freqdata1.label(3);

gdata.g1_multi_reg_conditional = ft_connectivityanalysis(grangercfg, freqdata1);
gdata.g2_multi_reg_conditional = ft_connectivityanalysis(grangercfg, freqdata2);

%% # Evaluation
% The label combinations are 6x2 cell arrays, containing all 2-permutations
% tuplets from the channels. How to interpret this?
% Is the combination a, b representing F<sub>a&rarr;b|c</sub>?
% Let's check this. In scenario 2, we should clearly see a higher causality
% from 1&rarr;3 | 2 than in scenario 1 of the differentially delayed
% drive. This corresponds to row 4 in the
% gdata.g1_multi_reg_conditional.labelcmb.
% So, let's compare the labelcmb 1, 3 in both scenarios:
%
scenario1_mean = mean(gdata.g1_multi_reg_conditional.grangerspctrm(4, :));
scenario2_mean = mean(gdata.g2_multi_reg_conditional.grangerspctrm(4, :));
