function test_bug811

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_checkdata spm2fieldtrip
% DATA private

ft_hastoolbox('spm8', 1);

filename = dccnpath('/project/3031000.02/test/bug811/cespm8_LD031286_DOELLER_20100223_01.mat');

D = spm_eeg_load(filename);
fnamedat = dccnpath('/project/3031000.02/test/bug811/cespm8_LD031286_DOELLER_20100223_01.dat');
E = clone(D, fnamedat);

data = spm2fieldtrip(E);
ft_checkdata(data, 'datatype', 'raw');
