function test_ft_artifact_xxx

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_artifact_clip ft_artifact_eog ft_artifact_muscle ft_artifact_zvalue ft_artifact_ecg ft_artifact_jump ft_artifact_threshold

% this script only checks whether the code runs without an error 
% it does not check whether the output is as expected

begsample = 1:600:6000;
begsample = begsample + 10000;
endsample = begsample + 600 - 1;
offset = zeros(size(begsample));

%%

cfg = [];
cfg.dataset     = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg.continuous  = 'yes';
cfg.trl         = [begsample(:) endsample(:) offset(:)];

% read the data into memory, this will be used later
data = ft_preprocessing(cfg);

cfg.artfctdef.ecg.channel        = 'MLT33'; % only one channel supported
cfg.artfctdef.ecg.feedback       = 'no';
cfg.artfctdef.clip.channel       = 'MLO';
cfg.artfctdef.threshold.channel  = 'ML';
cfg.artfctdef.zvalue.channel     = 'MLO';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.interactive = 'no';
cfg.artfctdef.eog.channel        = 'MLF';
cfg.artfctdef.eog.interactive    = 'no';
cfg.artfctdef.muscle.channel     = 'MLT';
cfg.artfctdef.muscle.interactive = 'no';
cfg.artfctdef.jump.channel       = 'MLO';
cfg.artfctdef.jump.interactive   = 'no';

[cfg, artifact1] = ft_artifact_ecg(cfg);
[cfg, artifact2] = ft_artifact_clip(cfg);
[cfg, artifact3] = ft_artifact_threshold(cfg);
[cfg, artifact4] = ft_artifact_zvalue(cfg);
[cfg, artifact5] = ft_artifact_eog(cfg);
[cfg, artifact6] = ft_artifact_muscle(cfg);
[cfg, artifact7] = ft_artifact_jump(cfg);

%%

cfg = [];
cfg.trl = [begsample(:) endsample(:) offset(:)];

cfg.artfctdef.ecg.channel        = 'MLT33'; % only one channel supported
cfg.artfctdef.ecg.feedback       = 'no';
cfg.artfctdef.clip.channel       = 'MLO';
cfg.artfctdef.threshold.channel  = 'ML';
cfg.artfctdef.zvalue.channel     = 'MLO';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.interactive = 'no';
cfg.artfctdef.eog.channel        = 'MLF';
cfg.artfctdef.eog.interactive    = 'no';
cfg.artfctdef.muscle.channel     = 'MLT';
cfg.artfctdef.muscle.interactive = 'no';
cfg.artfctdef.jump.channel       = 'MLO';
cfg.artfctdef.jump.interactive   = 'no';

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

