function test_bug1249

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_componentanalysis ft_rejectcomponent
% DATA private

load(dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf275.mat'));

cfg = [];
cfg.method = 'fastica';
cfg.numcomponent = 20;
cfg.channel = 'MEG';
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.component = 1:5;
datanew = ft_rejectcomponent(cfg, comp);

if size(datanew.grad.chanpos,1)~=numel(datanew.grad.label)
  error('labels in the grad structure are inconsistent with the chanpos');
end
