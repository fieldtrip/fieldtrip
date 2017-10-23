function test_bug46

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_connectivityanalysis univariate2bivariate ft_checkdata

% the problem reported was that ft_connectivityanalysis crashes with input is freq-data containing fourierspectra and two channels
% crash occurs in subfunction univariate2bivariate, because checkdata is called with cmbrepresentation 'sparse'
% whereas probably it should be 'cmbrepresentation is 'full', or a cmbindx should be given.
% error was reported for coh and psi as a method. unknown whether it also occurs for other metrics

% try to reproduce
rand('twister', 20070425);

freq       = [];
freq.label = {'a';'b'};
freq.freq  = 1:10;
freq.fourierspctrm = randn(9,2,10) + randn(9,2,10).*i;
freq.dimord = 'rpttap_chan_freq';
freq.cfg   = [];
freq.cumtapcnt = [3;3;3];

cfg        = [];
cfg.method = 'coh';
coh        = ft_connectivityanalysis(cfg, freq);

%this indeed leads to an error:
%
%??? Subscript indices must either be real positive integers or logicals.
%
%Error in ==> checkdata>fixcsd at 888
%    tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
%
%Error in ==> checkdata at 577
%    data = fixcsd(data, cmbrepresentation, channelcmb);
%
%Error in ==> ft_connectivityanalysis>univariate2bivariate at 992
%      data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb',
%      cmb);
%
%Error in ==> ft_connectivityanalysis at 138
%        [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm',
%        'crsspctrm', dtype, 0, cfg.channelcmb);
%

%uncommenting lines 124 to 126 seems to do the trick.
%I don't remember having commented out that part.

cfg         = [];
cfg.method  = 'coh';
cfg.complex = 'abs';
coh         = ft_connectivityanalysis(cfg, freq);

cfg         = [];
cfg.method  = 'psi';
cfg.bandwidth = 4;
psi         = ft_connectivityanalysis(cfg, freq);
