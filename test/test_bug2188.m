function test_bug2188

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY
% DATA private

%% Test read Localite markers

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2188/EEG_Markers.xml');
sens = ft_read_sens(filename);
