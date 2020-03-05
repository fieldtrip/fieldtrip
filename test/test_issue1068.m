function test_issue1068

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_checkdata ft_datatype makessense getdimord getdimsiz

%%
% this worked

freq= [];
freq.dimord = 'rpt_chan_freq';
freq.powspctrm = rand(10,3,1); % singleton dimension at the end
freq.label = {'1', '2', '3'};
freq.freq = 1;
ft_datatype(freq)
ft_checkdata(freq)


%%
% this failed

tlck = [];
tlck.dimord = 'rpt_chan_time';
tlck.trial = rand(10,3,1); % singleton dimension at the end
tlck.label = {'1' '2', '3'};
tlck.time = 1;
ft_datatype(tlck)
ft_checkdata(tlck)
