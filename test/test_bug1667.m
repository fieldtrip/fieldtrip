function test_bug1667

% WALLTIME 00:20:00
% MEM 4gb

% TEST ft_read_data
% TEST ft_fetch_data

% Improvements:
% 1-make buffering from multiple sources possible
% 2-if part of the trial is already buffered, only fetch what is needed

datainfo = ref_datasets;

datanr = 13;
datafile1 = fullfile(datainfo(datanr).origdir,'original',datainfo(datanr).type,datainfo(datanr).datatype,'/',datainfo(datanr).filename);
datanr = datanr+1;
datafile2 = fullfile(datainfo(datanr).origdir,'original',datainfo(datanr).type,datainfo(datanr).datatype,'/',datainfo(datanr).filename);
hdr1 = ft_read_header(datafile1);
hdr2 = ft_read_header(datafile2);
% fetch something
tic;
dat1 = ft_read_data(datafile1, 'cache', true, 'begsample', 1, 'endsample', 1000, 'checkboundary', false, 'header', hdr1);
t11 = toc;
tic;
dat2 = ft_read_data(datafile2, 'cache', true, 'begsample', 1, 'endsample', 1000, 'checkboundary', false, 'header', hdr1);
t12 = toc;
% fetch something that is partly overlapping
tic;
dat1 = ft_read_data(datafile1, 'cache', true, 'begsample', 500, 'endsample', 1500, 'checkboundary', false, 'header', hdr1);
t21 = toc;
tic;
dat2 = ft_read_data(datafile2, 'cache', true, 'begsample', 500, 'endsample', 1500, 'checkboundary', false, 'header', hdr1);
t22 = toc;
% fetch something completely else
tic;
dat1 = ft_read_data(datafile1, 'cache', true, 'begsample', 2000, 'endsample', 3000, 'checkboundary', false, 'header', hdr1);
t31 = toc;
tic;
dat2 = ft_read_data(datafile2, 'cache', true, 'begsample', 2000, 'endsample', 3000, 'checkboundary', false, 'header', hdr1);
t32 = toc;

% compare timing, the latter should be way slower, if not ERROR!!!11
