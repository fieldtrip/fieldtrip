function test_bug1746

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_sourceanalysis test_bug1746 ft_prepare_leadfield

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_singleshell.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/freq/meg/freq_mtmfft_powandcsd_ctf275'));
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/source/meg/source_grid_mtmfft_fourier_trl_DICS_fixedori_ctf275'));

sourcemodel = ft_source2grid(source);
sourcemodel.inside = sourcemodel.inside(1:10);
sourcemodel.outside = setdiff(1:prod(sourcemodel.dim),sourcemodel.inside);
meglabel    = ft_channelselection('MEG', freq.label);

cfg           = [];
cfg.headmodel = vol;
cfg.grid      = sourcemodel;
cfg.channel   = meglabel;
leadfield1    = ft_prepare_leadfield(cfg, freq);
assert(isfield(leadfield1, 'label'));

shuffle     = randperm(numel(meglabel));
cfg.channel = meglabel(shuffle);
leadfield2  = ft_prepare_leadfield(cfg, freq);
assert(norm(leadfield1.leadfield{4}(:)-leadfield2.leadfield{4}(:))./norm(leadfield1.leadfield{4})<10*eps);

% OBSERVATION: giving a different order of channels to
% ft_prepare_leadfield does not lead to a different ordering of the
% channels in the output. This is something which I don't worry about now
% because the leadfield.label is consistent with the data in
% leadfield.leadfield.

cfg = [];
cfg.method      = 'dics';
cfg.headmodel   = vol;
cfg.grid        = leadfield1;
cfg.frequency   = 5;
cfg.dics.lambda = '10%';
cfg.channel     = meglabel;
source1 = ft_sourceanalysis(cfg, freq);

cfg.channel = meglabel(1:2:end);
source2 = ft_sourceanalysis(cfg, freq);
