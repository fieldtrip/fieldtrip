clear all;

load('tmp_data_small.mat', 'data');
load('tmp_cfgart.mat', 'cfgart');

% Note that I marked artifacts in trials 1, 3, 5, 7
% 1: from the start of the trial to mid-way
% 3: in the middle of the trial
% 5: starting mid-way till the end of the trial
% 7: a complete trial
ft_databrowser(cfgart, data);

% Now reject those artifacts by filling them with NaNs
cfg = cfgart;
cfg.artfctdef.reject = 'nan';
dataClean = ft_rejectartifact(cfg, data);

% And inspect the result
ft_databrowser([], dataClean);

% The result is that trials 1, 3 and 5 are treated the way I would expect.
% Trial 7 however is still entirely present; it's neither removed, nor
% filled with NaNs.


