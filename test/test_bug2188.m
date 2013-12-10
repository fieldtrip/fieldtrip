function test_bug2188

% MEM 1500mb
% WALLTIME 00:03:02

%% Test read Localite markers

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2188/EEG_Markers.xml');
sens = ft_read_sens(filename);
