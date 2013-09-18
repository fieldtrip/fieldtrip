function test_bug2188

%% Test read Localite markers

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2188/EEG_Markers.xml');
sens = ft_read_sens(filename);