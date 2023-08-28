function test_issue1972

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_definetrial ft_read_event read_trigger
% DATA private

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/issue1972.fif');
values   = [14593 14594 14600];

% The reported problem by Mikkel, is that the above dataset does not
% return the correct trl matrix, unless the events are explicitly read with
% a 'detectflank', 'up' instruction. The fif-file contains both ST101 and
% ST102 trigger channels. Note that not all ST101/102 containing datasets
% might have this problem, so there might be an interaction between data
% and code. The working hypothesis is that the automatic detectflank
% specification in read_trigger is not sufficiently robust.

eve = ft_read_event(filename, 'detectflank', 'up');
           
cfg = [];
cfg.event               = eve;
cfg.dataset             = filename;
cfg.trialfun            = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'STI101';
cfg.trialdef.prestim    = 2.5;
cfg.trialdef.poststim   = 2.5;
cfg.trialdef.eventvalue = values;

cfg1 = ft_definetrial(cfg);

cfg = [];
% cfg.event               = eve;
cfg.dataset             = filename;
cfg.trialfun            = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'STI101';%'Trigger'
cfg.trialdef.prestim    = 2.5;
cfg.trialdef.poststim   = 2.5;
cfg.trialdef.eventvalue = values;%-14592;

cfg2 = ft_definetrial(cfg);

cfg.trialdef.detectflank = 'up';

cfg3 = ft_definetrial(cfg);

assert(isequal(cfg1.trl, cfg3.trl));
assert(isequal(cfg1.trl, cfg2.trl));
