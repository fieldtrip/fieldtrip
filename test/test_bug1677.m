function test_bug1677

% MEM 1500mb
% WALLTIME 00:10:00


% code contributed by Jasper Poort: thanks for that!

% problem with using a trial selection to do freqdescriptives on a subset of trials:

% example data
cfg             = [];
cfg.ntrials     = 10;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8    0    0 ; 
                        0  0.9  0.5 ;
                      0.4    0  0.5];                      
cfg.params(:,:,2) = [-0.5    0    0 ; 
                        0 -0.8    0 ; 
                        0    0 -0.2];                        
cfg.noisecov      = [ 0.3    0    0 ;
                        0    1    0 ;
                        0    0  0.2];
data              = ft_connectivitysimulation(cfg);

% freqanalysis tfr
cfg = [];
cfg.output      = 'fourier';
cfg.method      = 'mtmconvol';
cfg.keeptrials  = 'yes';
timunit         = 0.2;
raleigh         = 1/timunit;
cfg.foi         = [6*raleigh : raleigh : 20*raleigh];
numfoi          = length(cfg.foi);
cfg.t_ftimwin   = ones(1,numfoi) .* timunit;
cfg.tapsmofrq   = ones(1,numfoi) .*2 .* raleigh;
cfg.pad         = 10 .* timunit;
cfg.toi         = [0.2: timunit ./ 10 : 0.5];

freq = ft_freqanalysis(cfg, data);

% freqdescriptives
tmpcfg = [];
tmpcfg.jackknife = 'no';
tmpcfg.keeptrials = 'no';
%tmpcfg.channel    = freq.label([chnix1,chnix2])
tmpcfg.trials     = 2:5;
tmp = ft_freqdescriptives(tmpcfg,freq);

% error was reproduced and fixed by changing seloverdim
% no explicit assertion is made in this test function


%--------------------------------------------------------------------------
% The original error is pasted below

% Error using  + 
% Array dimensions must match for binary array op.
% 
% Error in ft_checkdata>fixcsd (line 740)
%       powspctrm = powspctrm + abs(data.fourierspctrm(p:ntap:end,:,:,:,:)).^2;
% 
% Error in ft_checkdata (line 646)
%     data = fixcsd(data, cmbrepresentation, channelcmb);
% 
% Error in ft_freqdescriptives (line 131)
% freq = ft_checkdata(freq, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
%  
% 740       powspctrm = powspctrm + abs(data.fourierspctrm(p:ntap:end,:,:,:,:)).^2

% the error seems to happen because the cumtapcnt seems to be updated after
% calling ft_selectdata in ft_freqdescriptives such that it no longer
% contains a row for every trial which causes nrpt =
% size(data.cumtapcnt,1); in ft_checkdata to return 1 instead of the number
% of trials.     
