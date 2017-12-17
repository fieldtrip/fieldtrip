function test_bug2040

% MEM 1500mb
% WALLTIME 00:10:00

% TEST nansum

[ftver, ftpath] = ft_version;

p  = tempname;
p1 = fullfile(p, 'mexfile');
p2 = fullfile(p, 'mfile');
p3 = fullfile(p, 'stats');

mkdir(p);
mkdir(fullfile(p1));
mkdir(fullfile(p2));
mkdir(fullfile(p3));

% separate the three implementations
cd(p1);
copyfile(fullfile(ftpath, 'src', ['nansum.' mexext]), '.');
cd(p2);
copyfile(fullfile(ftpath, 'src', 'nansum.m'), '.');
cd(p3);
% don't copy anything here, assume that the Mathworks signal processing toolbox is on the path

% make some data
n = 10;
x{1} = randn(n,1);
x{2} = x{1}';
x{3} = randn(n,n);
x{4} = randn(n,n,n);
x{5} = randn(n,n,n,n);

originalpath = path;
restoredefaultpath

% compute the solution for the three implementations
cd(p1);
s1 = cellfun(@nansum, x, 'UniformOutput', false);
cd(p2);
s2 = cellfun(@nansum, x, 'UniformOutput', false);
cd(p3);
s3 = cellfun(@nansum, x, 'UniformOutput', false);

% the solutions should all be equal
assert(isequal(s1, s2));
assert(isequal(s2, s3));

path(originalpath);

