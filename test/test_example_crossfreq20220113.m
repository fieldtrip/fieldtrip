function test_example_crossfreq

% MEM 4gb
% WALLTIME 00:10:00

%
%% Cross-frequency analysis
%
% There are several ways in which cross-frequency interactions might occur. In [Jensen and Colgin TICS 2007](http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6VH9-4NWNF64-3&_user=668715&_coverDate=07%2F31%2F2007&_rdoc=1&_fmt=&_orig=search&_sort=d&view=c&_acct=C000036278&_version=1&_urlVersion=0&_userid=668715&md5=2946af6effbd30ccf5d896ddfa6fa75c) different principles of cross-frequency interactions are shown in Figure 1.
%
%
% With the **[ft_freqsimulation](https://github.com/fieldtrip/fieldtrip/blob/release/ft_freqsimulation.m)** function you can generate simulated data in FieldTrip format which the different types of cross-frequency interactions. The different methods are:
%
%* [phalow_amphigh (is phase to power in Jensen and Colgin)](/example/crossfreq/phalow_amphigh)
%* [amplow_amphigh (is power to power in Jensen and Colgin](/example/crossfreq/amplow_amphigh)
%* [phalow_freqhigh (is phase to frequency in Jensen and Colgin)](/example/crossfreq/phalow_freqhigh)
