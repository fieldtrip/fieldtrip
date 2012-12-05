function test_bug1870

% TEST test_bug1870
% TEST test_ft_megplanar ft_megplanar ft_datatype_sens ft_compute_leadfield

% this test is basically a small section of test_ft_megplanar
% the input data is consistent, but along the way the grad structure gets screwed up

load /home/common/matlab/fieldtrip/data/test/bug1870.mat
data5 = ft_megplanar(cfg, data);
