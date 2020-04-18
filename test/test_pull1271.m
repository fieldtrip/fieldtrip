function test_pull1271

% MEM 5gb
% WALLTIME 00:20:00
% DEPENDENCY ft_read_data ft_read_header ft_read_event

% the DCCNPATH function will take the file from the present working directory if possible
% otherwise it will search on /home/common/matlab (a DCCN specific linux network drive)
% and if that fails it will search on H:\common\matlab (the windows version of that network drive)
%
% cd '/Users/alan/Dropbox/Lab-BML/Code/[2020.04]Fieldtrip_Neuroomega_nsx_PR';

fname_nsx         = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/blackrock/test_data_blackrock.ns2');
fname_neuroomega  = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/neuroomega/LT1D0.000F0000.mat');

%% loading blackrock header and data
hdr_nsx=ft_read_header(fname_nsx);
ft_read_data(fname_nsx,'header',hdr_nsx);

%% loading neuroomega header, events and data
hdr_neuroomega=ft_read_header(fname_neuroomega,'chantype','micro');
ft_read_data(fname_neuroomega,'header',hdr_neuroomega);
ft_read_event(fname_neuroomega);
