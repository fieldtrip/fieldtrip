function test_bug1807

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_senstype bti2grad ft_datatype_sens ft_read_header ft_read_sens

dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/bti248grad/e,rfhp1.0Hz,COH');

hdr  = ft_read_header(dataset);
grad = ft_read_sens(dataset);

assert(ft_senstype(hdr, 'bti248grad'));
assert(ft_senstype(grad, 'bti248grad'));

