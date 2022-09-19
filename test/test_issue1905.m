function test_issue1905

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY data2bids

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1905'));

% this contains the filename, and cfg.electrodes
load cfg

% specify the location for the BIDS output
cfg.bidsroot = tempname();

data2bids(cfg);

% clean up the BIDS output
rmdir(cfg.bidsroot, 's');
