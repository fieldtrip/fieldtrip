function test_pull2043

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY data2bids

global ft_default
ft_default.dccnpath = '/Volumes/SamsungT7/fieldtrip';

ctfdataset = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2043/emptyroom_Noise_20220613_01.ds');
neuromagdataset = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2043/empty_room_before_tsss_corr95.fif');

%%

% bidsroot = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/pull2043'), 'bids');
bidsroot = tempdir;

%%

cfg = [];
cfg.bidsroot = bidsroot;
cfg.dataset = ctfdataset;
cfg.datatype = 'meg';
cfg.method = 'copy';
cfg.sub = 'ctf';
cfg.ses = 'meg';
cfg.task = 'emptyroom';

cfg.participants.age = nan;
cfg.participants.sex = nan;

cfg.scans.acq_time = datetime('yesterday');

data2bids(cfg)

%%

gradfile = fullfile(bidsroot, 'sub-ctf', 'ses-meg', 'meg', 'sub-ctf_ses-meg_task-emptyroom_meg.ds');

grad = ft_read_sens(gradfile, 'coordsys', 'head', 'coilaccuracy', []);
assert(strcmp(grad.coordsys, 'unknown')) % "Other" from the coordsystem_json becomes "unknown"
assert(strcmp(grad.unit, 'cm')) % the cm units from the coordsystem_json should be enforced

grad = ft_read_sens(gradfile, 'coordsys', 'head', 'coilaccuracy', 2);
assert(strcmp(grad.coordsys, 'unknown')) % "Other" from the coordsystem_json becomes "unknown"
assert(strcmp(grad.unit, 'cm')) % the cm units from the coordsystem_json should be enforced

grad = ft_read_sens(gradfile, 'coordsys', 'dewar');
assert(strcmp(grad.coordsys, 'dewar')) % this returns the sensors in dewar coordinates, regardless of the coordsystem_json
assert(strcmp(grad.unit, 'cm')) % the cm units from the coordsystem_json should be enforced


%%

cfg = [];
cfg.bidsroot = bidsroot;
cfg.dataset = neuromagdataset;
cfg.datatype = 'meg';
cfg.method = 'copy';
cfg.sub = 'neuromag';
cfg.ses = 'meg';
cfg.task = 'emptyroom';

cfg.participants.age = nan;
cfg.participants.sex = nan;

cfg.scans.acq_time = datetime('today');

data2bids(cfg)

%%

gradfile = fullfile(bidsroot, 'sub-neuromag', 'ses-meg', 'meg', 'sub-neuromag_ses-meg_task-emptyroom_meg.fif');

grad = ft_read_sens(gradfile, 'coordsys', 'head', 'coilaccuracy', []);
assert(strcmp(grad.coordsys, 'unknown')) % "Other" from the coordsystem_json becomes "unknown"
assert(strcmp(grad.unit, 'm')) % the m units from the coordsystem_json should be enforced

grad = ft_read_sens(gradfile, 'coordsys', 'head', 'coilaccuracy', 2);
assert(strcmp(grad.coordsys, 'unknown')) % "Other" from the coordsystem_json becomes "unknown"
assert(strcmp(grad.unit, 'm')) % the m units from the coordsystem_json should be enforced

grad = ft_read_sens(gradfile, 'coordsys', 'dewar');
assert(strcmp(grad.coordsys, 'dewar')) % this returns the sensors in dewar coordinates, independent of the coordsystem_json
assert(strcmp(grad.unit, 'm')) % the m units from the coordsystem_json should be enforced

%%

rmdir(bidsroot, 's');
