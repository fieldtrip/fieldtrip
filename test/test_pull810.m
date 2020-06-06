function pull810

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_fetch_data ft_databrowser ft_appenddata

% Contact:
% Johannes Algermissen
% j.algermissen@donders.ru.nl
% project 3017042.02
% PI Hanneke den Ouden

% contains 6 blocks to be concatenated
load(dccnpath('/home/common/matlab/fieldtrip/data/test/pull810.mat'))

%%
% this makes it easier to debug
for i=1:6
  for j=1:numel(blockdata{i}.trial)
    blockdata{i}.trial{j}(:) = i;
  end
end

%% Concatenate blocks:
% a) Link normally:
cfg = [];
data = ft_appenddata(cfg, blockdata{1:6}); % concatenate blocks
% --> when plotted: gives weird first trial

%%
% b) Via looping:
data      = blockdata{1};
for iBlock = 2:6
  cfg = [];
  data = ft_appenddata(cfg, data, blockdata{iBlock}); % concatenate blocks
end
% --> same as a): when plotted: gives weird first trial

%%
try
  cfg = [];
  cfg.channel = 1;
  cfg.continuous = 'no';
  ft_databrowser(cfg, data);
  failed = false;
catch
  failed = true;
end
% the default is cfg.allowoverlap='no'
assert(failed, 'this should have resulted in an error');

cfg = [];
cfg.channel = 1;
cfg.continuous = 'no';
cfg.allowoverlap = 'yes';
ft_databrowser(cfg, data);

%%
% using cfg.keepsampleinfo option
data = blockdata{1};
for iBlock = 2:6
  cfg = [];
  cfg.keepsampleinfo = 'no';
  data = ft_appenddata(cfg, data, blockdata{iBlock}); % concatenate blocks
end

%%

cfg = [];
cfg.channel = 1;
cfg.continuous = 'no';
ft_databrowser(cfg, data);
