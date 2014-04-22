function test_bug2188

% MEM 1500mb
% WALLTIME 00:10:00

%% Test read Localite markers

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2188/EEG_Markers.xml');
sens = ft_read_sens(filename);
