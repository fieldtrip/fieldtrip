function test_ft_conjunctionanalysis

% MEM 1500mb
% WALLTIME 00:10:00

% FT_CONJUNCTIONANALYSIS finds the minimum statistic common across two or 
% more contrasts, i.e. data following ft_xxxstatistics. Furthermore, it 
% finds the overlap of sensors/voxels that show statistically significant 
% results (a logical AND on the mask fields).
%
% Alternatively, it finds minimalistic mean power values in the
% input datasets. Here, a type 'relative change' baselinecorrection
% prior to conjunction is advised.
% A. Stolk

% chan-time data
timelock1.label = {'chan1';'chan2'};
timelock1.time  = 1:10;
timelock1.dimord = 'chan_time';
timelock1.stat  = randn(2,10);
timelock2 = timelock1;
timelock2.stat  = randn(2,10);

stat = ft_conjunctionanalysis([], timelock1, timelock2);

timelock3 = timelock2;
timelock3.stat  = randn(2,10);

stat = ft_conjunctionanalysis([], timelock1, timelock2, timelock3);

% 2 sinuses - for fun
timelock1.label = {'chan1';'chan2'};
timelock1.time  = [1/1000:1/1000:1];
timelock1.dimord = 'chan_time';
timelock1.stat  = sin(2*pi*2*timelock1.time);
timelock2 = timelock1;
timelock2.stat  = sin(2*pi*2*timelock1.time+200);
stat = ft_conjunctionanalysis([], timelock1, timelock2);

% chan-freq data
freq1.label = {'chan1';'chan2'};
freq1.freq  = 1:10;
freq1.dimord = 'chan_freq';
freq1.stat  = randn(2,10);
freq2 = freq1;
freq2.stat  = randn(2,10);

stat = ft_conjunctionanalysis([], freq1, freq2);

freq3 = freq2;
freq3.stat  = randn(2,10);

stat = ft_conjunctionanalysis([], freq1, freq2, freq3);

% create some chan-freq-time data
freq1.label = {'chan1';'chan2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.stat  = randn(2,10,5);
freq2 = freq1;
freq2.stat  = randn(2,10,5);

stat = ft_conjunctionanalysis([], freq1, freq2);

freq3 = freq2;
freq3.stat  = randn(2,10,5);

stat = ft_conjunctionanalysis([], freq1, freq2, freq3);

% source data
%source1.dim = [50 1];
source1.inside = [1:25];
source1.outside = [26:50];
source1.pos = randn(50,3);
source1.stat = randn(50,1);
source1.prob = randn(50,1);
source1.mask = ones(50,1);
source2 = source1;
source2.stat = randn(50,1);
source2.prob = randn(50,1);

stat = ft_conjunctionanalysis([], source1, source2);

source3 = source2;
source3.stat = randn(50,1);
source3.prob = randn(50,1);

stat = ft_conjunctionanalysis([], source1, source2, source3);
