function test_bug3077

% MEM 1gb
% WALLTIME 00:20:00

% TEST ft_rejectcomponent

% create some data
data = [];
data.trial = {randn(2,1000)};
data.time  = {(0:999)./1000};
data.label = {'chan01';'chan02'};
data.cfg = [];

filename1 = [tempname,'.mat'];
filename2 = [tempname,'.mat'];
filename3 = [tempname,'.mat'];

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.component = 2;
clean1 = ft_rejectcomponent(cfg, comp);
clean2 = ft_rejectcomponent(cfg, comp, data);

save(filename1, 'comp');
save(filename2, 'data');

cfg.inputfile = {filename1, filename2};
cfg.outputfile = filename3;
ft_rejectcomponent(cfg);

clean3 = load(filename3);
clean3 = clean3.data;

% there are small numerical differences
assert(norm(clean1.trial{1}-clean2.trial{1})<100*eps);
assert(norm(clean1.trial{1}-clean3.trial{1})<100*eps);
assert(norm(clean2.trial{1}-clean3.trial{1})<100*eps);

delete(filename1);
delete(filename2);
delete(filename3);

