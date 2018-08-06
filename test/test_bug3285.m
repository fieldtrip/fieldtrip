function test_bug3285

% WALLTIME 00:10:00
% MEM 2gb

% load('SubjectUCI29_data.mat', 'data');
cd(dccnpath('/home/common/matlab/fieldtrip/data/test'));
load bug3285.mat

%%

cfg               = [];
cfg.channel       = ft_channelselection({'LPG*', 'LTG*'}, data.label);
cfg.reref         = 'yes';
cfg.refchannel    = 'all';
cfg.updatesens    = 'no';
reref_grids = ft_preprocessing(cfg, data);

% CORRECT: chansel is not applied to elec struc
assert(isequal(data.elec, reref_grids.elec));


%%

depths = {'RAM*', 'RHH*', 'RTH*', 'ROC*', 'LAM*', 'LHH*', 'LTH*'};
for d = 1:numel(depths)
  cfg                    = [];
  cfg.channel            = ft_channelselection(depths{d}, data.label);
  cfg.montage.labelold   = cfg.channel;
  cfg.montage.labelnew   = strcat(cfg.channel(1:end-1),'-',cfg.channel(2:end));
  cfg.montage.tra        = [
    1    -1     0     0     0     0     0     0
    0     1    -1     0     0     0     0     0
    0     0     1    -1     0     0     0     0
    0     0     0     1    -1     0     0     0
    0     0     0     0     1    -1     0     0
    0     0     0     0     0     1    -1     0
    0     0     0     0     0     0     1    -1
    ];
  cfg.updatesens = 'yes';
  reref_depths{d} = ft_preprocessing(cfg, data);
  
  %?CORRECT: montage is applied to elec struc, i.e. tra is updated
  assert(~isequal(data.elec, reref_depths{d}.elec));
  
end



%%


cfg             = [];
% cfg.appendsens  = 'no'; this is the default
reref_all = ft_appenddata(cfg, reref_grids, reref_depths{:});

% CORRECT: input elec is not consistent, hence output shoudl not have an elec
assert(~isfield(reref_all, 'elec'));


%%

try
  cfg             = [];
  cfg.appendsens  = 'yes';
  reref_all = ft_appenddata(cfg, reref_grids, reref_depths{:});
  
  % FIXME: it is not clear what should happen here
end



