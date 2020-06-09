function test_ft_selectdata_savevar

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

%%

timelock1 = [];
timelock1.label = {'1' '2', '3'};
timelock1.time  = 1:5;
timelock1.dimord = 'chan_time';
timelock1.avg = randn(3,5);
timelock1.cfg = 'initial 1'; % note that this is a string rather than a struct

timelock2 = [];
timelock2.label = {'2' '3', '4'};
timelock2.time  = 2:6;
timelock2.dimord = 'chan_time';
timelock2.avg = randn(3,5);
timelock2.cfg = 'initial 2'; % note that this is a string rather than a struct

%%

cfg = [];
filename1 = [tempname,'.mat'];
filename2 = [tempname,'.mat'];
cfg.outputfile = {filename1, filename2};

[timelock1sel, timelock2sel] = ft_selectdata(cfg, timelock1, timelock2);

%%

% initial sanity check: it should return the intersection, i.e. 2 channels and 4 timepoints
assert(numel(timelock1sel.label)==2);
assert(numel(timelock2sel.label)==2);
assert(numel(timelock1sel.time)==4);
assert(numel(timelock2sel.time)==4);

% check the history postamble
assert(isfield(timelock1sel, 'cfg') && isstruct(timelock1sel.cfg));
assert(isfield(timelock2sel, 'cfg') && isstruct(timelock2sel.cfg));

% check the provenance postamble
assert(isfield(timelock1sel.cfg, 'previous') && isequal(timelock1sel.cfg.previous, 'initial 1'));
assert(isfield(timelock2sel.cfg, 'previous') && isequal(timelock2sel.cfg.previous, 'initial 2'));

% check the savevar postamble
assert(exist(filename1, 'file')~=0);
assert(exist(filename2, 'file')~=0);
tmp1 = load(filename1);
tmp2 = load(filename2);
delete(filename1);
delete(filename2);
assert(isequal(tmp1.data, timelock1sel));
assert(isequal(tmp2.data, timelock2sel));
