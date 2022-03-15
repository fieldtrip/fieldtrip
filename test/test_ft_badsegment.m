function test_ft_badsegment

% WALLTIME 00:10:00
% MEM 2gb

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');

ft_debug off
interactive = false;

%%

cfg = [];
cfg.dataset = dataset;
cfg.channel = 'MEG';
cfg.trl = [1 7*3*300 0]; % first 7 trials, each of 3 seconds
cfg.continuous = 'yes';
cfg.demean = 'yes';
datacnt = ft_preprocessing(cfg);

cfg = [];
cfg.length = 1;
datatrl = ft_redefinetrial(cfg, datacnt);

%%

cfg = [];
cfg.method = 'template';
neighbours = ft_prepare_neighbours(cfg, datatrl);

%%

if interactive
  % this can be used to determine the thresholds
  cfg = [];
  cfg.method = 'summary';
  cfg.neighbours = neighbours;
  cfg.metric = 'var';
  data_rejected = ft_rejectvisual(cfg, datatrl);
end

%%
% each method has its own optimal threshold
% FT_REJECTVISUAL can be used to determine the appropriate value

cfg = [];
cfg.feedback = 'yes';
cfg.neighbours = neighbours;

cfg.metric = 'var';
cfg.threshold = 1.5e-25;

% cfg.metric = 'std';
% cfg.threshold = 4.2e-13;
%
% cfg.metric = 'max';
% cfg.threshold = 1.5e-12;
%
% cfg.metric = 'min';
% cfg.threshold = -12e-13;
%
% cfg.metric = 'range';
% cfg.threshold = 3000e-15;
%
% cfg.metric = 'zvalue';
% cfg.threshold = 3;
% 
% cfg.metric = 'neighbexpvar';
% cfg.threshold = 0.5;
%
% cfg.metric = 'neighbstdratio';
% cfg.threshold = 0.5;
%
% cfg.metric = 'neighbcorr';
% cfg.threshold = 0.6;
%
% these options can be varied for neighbstdratio and neighbcorr
% cfg.nbdetect = 'median';
% cfg.nbdetect = 'any';
% cfg.nbdetect = 'most';
% cfg.nbdetect = 'all';

outcfg = ft_badsegment(cfg, datatrl);

%%

cfg = [];
cfg.artfctdef.badsegment.artifact = outcfg.artfctdef.badsegment.artifact;
cfg.continuous = 'yes';
cfg.blocksize = 21;
cfg.colorgroups = ones(size(datatrl.label));
cfg.linecolor = 'br';
cfg.layout = 'CTF151.lay';
ft_databrowser(cfg, datatrl);

%%

cfg = [];
cfg.artfctdef.badsegment.artifact = outcfg.artfctdef.badsegment.artifact;
data_rejected = ft_rejectartifact(cfg, datatrl);
