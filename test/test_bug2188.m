function test_bug2188

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

%% Test read Localite markers

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2188/EEG_Markers.xml');
sens = ft_read_sens(filename);
