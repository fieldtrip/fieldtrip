function test_bug3219

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_rejectartifact
% TEST ft_rejectvisual

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3218.mat'), 'data'); % on purpose

% Note that I marked artifacts in trials 1, 3, 5, 7
% 1: from the start of the trial to mid-way
% 3: in the middle of the trial
% 5: starting mid-way till the end of the trial
% 7: a complete trial
% cfgart = ft_databrowser([], data);
cfgart.artfctdef.visual.artifact = [257371      258063
      261434      261922
      267309      268020
      269681      271240];

% Now reject those artifacts by filling them with NaNs
cfg = cfgart;
cfg.artfctdef.reject = 'nan';
dataClean = ft_rejectartifact(cfg, data);

% And inspect the result
% ft_databrowser([], dataClean);

hasnans = false(numel(dataClean.trial),1);
for k = 1:numel(dataClean.trial)
  % global check for nans in the trials
  hasnans(k,1) = any(isnan(dataClean.trial{k}(:)));
end

assert(isequal(hasnans', [true false true false true false true false]));

% The result is that trials 1, 3 and 5 are treated the way I would expect.
% Trial 7 however is still entirely present; it's neither removed, nor
% filled with NaNs.


