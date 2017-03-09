function test_bug1618

% MEM 1500mb
% WALLTIME 00:10:00


filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1618/bug1618.dat');

h = ft_read_header(filename);

disp(h)
disp(h.chanunit)
disp(h.chantype)

X = ft_read_data(filename);
assert(~any(isnan(X(:))), 'NaN values are not expected for this dataset!');
