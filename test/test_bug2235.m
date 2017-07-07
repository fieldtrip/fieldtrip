function test_bug2235

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_denoise_synthetic

fname = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2235');
load(fname);

cfg          = [];
cfg.gradient = 'G3BR';
da_data      = ft_denoise_synthetic(cfg, a_data);

% check label order in the data
[da,a]=match_str(da_data.label,a_data.label);
assert(all(da==a),'order of the labels is different in input and output data');
[da,a]=match_str(da_data.grad.label,a_data.grad.label);

% assert(all(da==a),'order of the labels is different in gradiometer description of input and output data');

% check whether the gradiometer arrays are internally consistent
[da,a] = match_str(da_data.grad.label,a_data.grad.label);
tra_a  = a_data.grad.tra(274:end,:);
tra_da = da_data.grad.tra(274:end,:);
assert(isequal(tra_a(a(274:end)-273,:),tra_da(da(274:end)-273,:)));
assert(isequal(da_data.grad.chanpos(da,:),a_data.grad.chanpos(a,:)));
assert(isequal(da_data.grad.chanori(da,:),a_data.grad.chanori(a,:)));
assert(isequal(da_data.grad.coilpos,a_data.grad.coilpos));
assert(isequal(da_data.grad.coilori,a_data.grad.coilori));

% and here is where it breaks
assert(isequal(da_data.grad.chantype(da),a_data.grad.chantype(a)));

