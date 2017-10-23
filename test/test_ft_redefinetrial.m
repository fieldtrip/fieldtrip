function test_ft_redefinetrial

% MEM 1500mb
% WALLTIME 00:10:00


load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151'));

data.trialinfo = (1:10)';

cfg        = [];
cfg.trials = 1:2:10;
data2      = ft_redefinetrial(cfg, data);
assert(all(data2.trialinfo==[1 3 5 7 9]'));

cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
data3        = ft_redefinetrial(cfg, data);
assert(all(data3.sampleinfo(:,1)==[51:300:2751]'));

cfg        = [];
cfg.toilim = [0.2 0.8];
data4      = ft_redefinetrial(cfg, data);
assert(all(data4.sampleinfo(:,1)==[61:300:2761]') && all(data4.sampleinfo(:,2)==[241:300:2941]'));

cfg        = [];
cfg.length = 0.2;
data5      = ft_redefinetrial(cfg, data);
assert(numel(data5.trial)==50 && all(data5.trialinfo(1:5)==1));

% this should give an error
% cfg        = [];
% cfg.trl    = [1 600 0];
% data6      = ft_redefinetrial(cfg, data);
% assert(~isfield(data6, 'trialinfo'));

% here the trialinfo should be added
data = rmfield(data, 'trialinfo');
cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
cfg.trl(:,4) = 1:10';
data7        = ft_redefinetrial(cfg, data);
assert(all(data7.trialinfo(:,1)==(1:10)'));

% here the user-specified trialinfo should be added
data.trialinfo = (1:10)';
cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
cfg.trl(:,4) = (11:20)';
data8        = ft_redefinetrial(cfg, data);
assert(all(data8.trialinfo(:,1)==(11:20)'));

% here the old trialinfo should be added
data.trialinfo = (1:10)';
cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
data8        = ft_redefinetrial(cfg, data);
assert(all(data8.trialinfo(:,1)==(1:10)'));

