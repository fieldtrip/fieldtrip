function test_ft_artifact_xxx

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_artifact_clip ft_artifact_eog ft_artifact_muscle ft_artifact_zvalue ft_artifact_ecg ft_artifact_jump ft_artifact_threshold

% this script only checks whether the code runs without an error 
% it does not check whether the output is as expected

begsample = 1:600:6000;
begsample = begsample + 10000;
endsample = begsample + 600 - 1;
offset = zeros(size(begsample));

cfg = [];
cfg.dataset     = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.continuous  = 'yes';
cfg.trl         = [begsample(:) endsample(:) offset(:)];

% read the data into memory, this will be used later
data = ft_preprocessing(cfg);

cfg.artfctdef.ecg.channel       = 'MLT33'; % only one channel supported
cfg.artfctdef.clip.channel      = 'MLO';
cfg.artfctdef.threshold.channel = 'ML';
cfg.artfctdef.zvalue.channel    = 'MLO';
cfg.artfctdef.eog.channel       = 'MLF';
cfg.artfctdef.muscle.channel    = 'MLT';
cfg.artfctdef.jump.channel      = 'MLO';
cfg.artfctdef.zvalue.cutoff     = 4;
cfg.artfctdef.ecg.feedback      = 'no';
cfg.artfctdef.muscle.feedback   = 'no';
cfg.artfctdef.zvalue.feedback   = 'no';
cfg.artfctdef.zvalue.feedback   = 'no';

[cfg, artifact] = ft_artifact_ecg(cfg);
[cfg, artifact] = ft_artifact_clip(cfg);
[cfg, artifact] = ft_artifact_threshold(cfg);
[cfg, artifact] = ft_artifact_zvalue(cfg);
[cfg, artifact] = ft_artifact_eog(cfg);
[cfg, artifact] = ft_artifact_muscle(cfg);
[cfg, artifact] = ft_artifact_jump(cfg);

cfg1 = cfg;
clear cfg;


cfg = [];
cfg.trl         = [begsample(:) endsample(:) offset(:)];
cfg.artfctdef.ecg.channel       = 'MLT33'; % only one channel supported
cfg.artfctdef.clip.channel      = 'MLO';
cfg.artfctdef.threshold.channel = 'ML';
cfg.artfctdef.zvalue.channel    = 'MLO';
cfg.artfctdef.eog.channel       = 'MLF';
cfg.artfctdef.muscle.channel    = 'MLT';
cfg.artfctdef.jump.channel      = 'MLO';
cfg.artfctdef.zvalue.cutoff     = 4;
cfg.artfctdef.ecg.feedback      = 'no';
cfg.artfctdef.muscle.feedback   = 'no';
cfg.artfctdef.zvalue.feedback   = 'no';
cfg.artfctdef.zvalue.feedback   = 'no';

cfg.artfctdef.eog.trlpadding = 0;
cfg.artfctdef.eog.artpadding = 0;
cfg.artfctdef.eog.fltpadding = 0;

cfg.artfctdef.muscle.trlpadding = 0;
cfg.artfctdef.muscle.artpadding = 0;
cfg.artfctdef.muscle.fltpadding = 0;

cfg.artfctdef.jump.trlpadding = 0;
cfg.artfctdef.jump.artpadding = 0;
cfg.artfctdef.jump.fltpadding = 0;

[cfg, artifact1] = ft_artifact_ecg(cfg, data);
[cfg, artifact2] = ft_artifact_clip(cfg, data);
[cfg, artifact3] = ft_artifact_threshold(cfg, data);
[cfg, artifact4] = ft_artifact_zvalue(cfg, data);
[cfg, artifact5] = ft_artifact_eog(cfg, data);
[cfg, artifact6] = ft_artifact_muscle(cfg, data);
[cfg, artifact7] = ft_artifact_jump(cfg, data);

