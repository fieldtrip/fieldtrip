function test_bug3144

% WALLTIME 00:10:00
% MEM 2gb

% TEST test_bug3144
% TEST ft_checkdata ft_datatype_sens

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3144.mat'))

ft_checkdata(data);

ft_datatype_sens(data.grad);
