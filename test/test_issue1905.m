function test_issue1905

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY data2bids
% DATA private

cd(dccnpath('/project/3031000.02/test/issue1905'));

% this contains the filename, and cfg.electrodes
load cfg

% specify the location for the BIDS output
cfg.bidsroot = tempname();

data2bids(cfg);

% clean up the BIDS output
rmdir(cfg.bidsroot, 's');
