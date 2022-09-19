function opt = options(opt)
% Get/Set options associated with the bayesFactor toolbox
% This function is called by bf internal functions. 
% User can fine tune options by first calling this function without 
% input arguments, then modifying the struct that is returned and then
% calling this function again, with the modified struct as its input.
% EG. To use 20 parpool workers by default, use:
% opt= bf.options
% opt.nrWorkers = 20;
% bf.options(opt);
%
% Please note that changes to the options are permanent (i.e. they will be
% used again in the next Matlab session on the same machine).


if nargin<1
%% Default parameters for the Monte Carlo integration (see internal.mcIntegrate) and parallel execution
opt.minG = 0.0001;  % Smalles g value (don't include zero).
opt.maxG = 10000;   % Largest g value to intergare
opt.stepG = 0.05;   % Step size for the g integration.
opt.nrSamples = 10000;  % How many g samples to draw for the MC integration
opt.nDimsForMC = 1; % If there are this many dimensions, use MC integration
opt.verbose = true; % Show messages
% 1 and 2D integration work fine with standard
% integral.m but 3d is slow and higher not
% possible.
opt.nrWorkers = 0;  % Set to zero to use serial execution (i.e. disable parfor loops in the code)   
end

if ~ispref('bayesFactor') || nargin>0
    % First run : setup pref
    % Or opt specified - save as pref for future use.
    setpref('bayesFactor','options',opt);
else 
    opt= getpref('bayesFactor','options');
end
