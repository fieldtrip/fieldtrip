% FT_PREAMBLE_RANDOMSEED is a helper script that records, and optionally
% sets, the random seed/state of the random number generators in
% RAND/RANDN/RANDI.  
% It calls RANDOMSEED which deals with the MATLAB version dependencies.
% If cfg.randomseed does not exist, it is set to empty, which indicates
% default behavior.  
% Default behavior is to only record the current state, not set it (except
% for Matlab versions prior to 7.3 which set it to a new random state).

cfg.randomseed=ft_getopt(cfg, 'randomseed', []);
ftFuncRandomseed=randomseed(cfg.randomseed);