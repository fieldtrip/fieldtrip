function test_pull2111

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_read_data, ft_read_header
% DATA private

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2110/data_test_pull2111.mat');

hdr = ft_read_header(filename, 'headerformat', 'neuroomega_mat', 'chantype', 'ecog');
data = ft_read_data(filename, 'header', hdr, 'dataformat', 'neuroomega_mat', 'headerformat', 'neuroomega_mat');
