function test_issue1089

% MEM 1gb
% WALLTIME 10:00:00
% DEPENDENCY ft_selectdata


%%
% make it really simple
%
% timelock1 = [];
% timelock1.label = {'1'};
% timelock1.time = 1:10;
% timelock1.avg = randn(1,10);
%
% timelock2 = [];
% timelock2.label = {'2'};
% timelock2.time = 11:15;
% timelock2.avg = randn(1,5);
%
% cfg = [];
% cfg.select = 'union';
% [sel1, sel2] = ft_selectdata(cfg, timelock1, timelock2);
%


%%
% make two datasets

freqa = [];
freqa.cfg = struct();
freqa.dimord = 'rpt_chan_freq_time';
freqa.label = cellfun(@num2str, num2cell(1:10), 'UniformOutput', false);
freqa.freq = 1:10;
freqa.time = 1:10;
freqa.powspctrm = randn(30, 10, 10, 10);
freqa.cumtapcnt = ones(30, 1);

% with zero overlap
freqb = [];
freqb.cfg = struct();
freqb.dimord = 'rpt_chan_freq_time';
freqb.label = cellfun(@num2str, num2cell(11:15), 'UniformOutput', false);
freqb.freq = 11:15;
freqb.time = 11:15;
freqb.powspctrm = randn(30, 5, 5, 5);
freqb.cumtapcnt = ones(30, 1);

runtest(freqa, freqb)

% with perfect overlap
freqb = [];
freqb.cfg = struct();
freqb.dimord = 'rpt_chan_freq_time';
freqb.label = cellfun(@num2str, num2cell(1:10), 'UniformOutput', false);
freqb.freq = 1:10;
freqb.time = 1:10;
freqb.powspctrm = randn(30, 10, 10, 10);
freqb.cumtapcnt = ones(30, 1);

runtest(freqa, freqb)

% with some overlap
freqb = [];
freqb.cfg = struct();
freqb.dimord = 'rpt_chan_freq_time';
freqb.label = cellfun(@num2str, num2cell(3:7), 'UniformOutput', false);
freqb.freq = 3:7;
freqb.time = 3:7;
freqb.powspctrm = randn(30, 5, 5, 5);
freqb.cumtapcnt = ones(30, 1);

runtest(freqa, freqb)

% with full overlap in channels and frequencies, and partial in time
freqb = [];
freqb.cfg = struct();
freqb.dimord = 'rpt_chan_freq_time';
freqb.label = cellfun(@num2str, num2cell(1:10), 'UniformOutput', false);
freqb.freq = 1:10;
freqb.time = 6:14;
freqb.powspctrm = randn(30, 10, 10, 10);
freqb.cumtapcnt = ones(30, 1);

runtest(freqa, freqb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runtest(freqa, freqb)

%%
% do the selection and averaging in one go

cfg = [];
cfg.channel          = 'all';
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.parameter        = 'powspctrm';
cfg.select           = 'union';

cfg.avgoverrpt       = 'yes';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
[freqa1r,freqb1r] = ft_selectdata(cfg, freqa, freqb);

cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'yes';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
[freqa1c,freqb1c] = ft_selectdata(cfg, freqa, freqb);

cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'no';
[freqa1t,freqb1t] = ft_selectdata(cfg, freqa, freqb);

cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'yes';
[freqa1f,freqb1f] = ft_selectdata(cfg, freqa, freqb);

% this is selecting without averaging, needed for comparison below
cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
[freqa1,freqb1] = ft_selectdata(cfg, freqa, freqb);


%%
% do the averaging after the selection

cfg = [];
cfg.select           = 'intersect'; % it is now a perfect intersection

cfg.avgoverrpt       = 'yes';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
[freqa2r,freqb2r] = ft_selectdata(cfg, freqa1, freqb1);
assert(isequaln(freqa1r.powspctrm, freqa2r.powspctrm))
assert(isequaln(freqb1r.powspctrm, freqb2r.powspctrm))

cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'yes';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
[freqa2c,freqb2c] = ft_selectdata(cfg, freqa1, freqb1);
assert(isequaln(freqa1c.powspctrm, freqa2c.powspctrm))
assert(isequaln(freqb1c.powspctrm, freqb2c.powspctrm))

cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'yes';
cfg.avgoverfreq      = 'no';
[freqa2t,freqb2t] = ft_selectdata(cfg, freqa1, freqb1);
assert(isequaln(freqa1t.powspctrm, freqa2t.powspctrm))
assert(isequaln(freqb1t.powspctrm, freqb2t.powspctrm))

cfg.avgoverrpt       = 'no';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'yes';
[freqa2f,freqb2f] = ft_selectdata(cfg, freqa1, freqb1);
assert(isequaln(freqa1f.powspctrm, freqa2f.powspctrm))
assert(isequaln(freqb1f.powspctrm, freqb2f.powspctrm))

