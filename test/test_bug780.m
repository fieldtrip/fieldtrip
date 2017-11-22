function test_bug780

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preproc_lowpassfilter ft_preproc_highpassfilter ft_preproc_bandpassfilter ft_preproc_bandstopfilter filter_with_correction

dat = randn(1,1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_preproc_highpassfilter(dat, 1000, 0.1, 9, [], [], 'split');
ft_preproc_highpassfilter(dat, 1000, 0.1, 9, [], [], 'reduce');
try
  ft_preproc_highpassfilter(dat, 1000, 0.1, 9, [], [], 'none');
  fail = false;
catch
  fail = true;
end
if ~fail
  error('this test should have failed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_preproc_lowpassfilter(dat, 1000, 0.1, 9, [], [], 'split');
ft_preproc_lowpassfilter(dat, 1000, 0.1, 9, [], [], 'reduce');
try
  ft_preproc_lowpassfilter(dat, 1000, 0.1, 9, [], [], 'none');
  fail = false;
catch
  fail = true;
end
if ~fail
  error('this test should have failed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_preproc_bandpassfilter(dat, 1000, [49 51], 9, [], [], 'split');
ft_preproc_bandpassfilter(dat, 1000, [49 51], 9, [], [], 'reduce');
try
  ft_preproc_bandpassfilter(dat, 1000, [49 51], 9, [], [], 'none');
  fail = false;
catch
  fail = true;
end
if ~fail
  error('this test should have failed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_preproc_bandstopfilter(dat, 1000, [49 51], 9, [], [], 'split');
ft_preproc_bandstopfilter(dat, 1000, [49 51], 9, [], [], 'reduce');
try
  ft_preproc_bandstopfilter(dat, 1000, [49 51], 9, [], [], 'none');
  fail = false;
catch
  fail = true;
end
if ~fail
  error('this test should have failed');
end
