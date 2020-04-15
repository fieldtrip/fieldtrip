% FT_EXAMPLE_READ_MEF3P0 an example to read MEF 2.1 data set into FieldTrip

%% clear all
clear

% look for session path of data
p = mfilename('fullpath');
sesspath = fullfile(fileparts(p), 'mef_2p1');

% set the password
password = struct('Subject', 'erlichda', 'Session', 'sieve', 'Data', '');

% get the object to deal with FieldTrip
mef_ft = MEFFieldTrip_2p1(sesspath, password);

%% now let's import two channels D_1 and F_1
% for the first 10 seconds
mef_ft.SelectedChannel = ["D_1" "F_1"];
mef_ft.StartEnd = [0 10];
mef_ft.SEUnit = 'Second';
[X, t] = mef_ft.importSession;
fs = mef_ft.SamplingFrequency;

figure
ph = plot(t/fs, X);
legend(ph, mef_ft.SelectedChannel)
xlim([0, 1])
xlabel('Time (s)')
ylabel('Amplitude')

%% let's demonstrate more interaction with FieldTrip. Let's import 10 seconds
% data at the beginning of the recording with an assumption that the
% trigger was at 0 second.
in_unit = 'second';
be_second = [0, 10, 0]; % 10-second time of data from the start
out_unit = 'index';
be_sample = mef_ft.SessionUnitConvert(be_second, in_unit, out_unit);

% read data header with ft_read_header
hdr = ft_read_header(sesspath, 'password', password);

% read data with ft_read_data()
dat = ft_read_data(sesspath, 'begsample', be_sample(1), 'endsample', be_sample(2),...
    'header', hdr, 'password', password);

% read data with ft_preprocessing
cfg = [];
cfg.dataset = sesspath;
cfg.password = password;
cfg.header = hdr;
cfg.trl = be_sample;
dat_ieeg = ft_preprocessing(cfg);

% plot the data
cfg.viewmode = 'vertical';
brwview = ft_databrowser(cfg, dat_ieeg);

%% Copyright 2020 Richard J. Cui
% Created: Thu 04/02/2020 11:30:27.749 AM
% Revision: 0.1  Date: Thu 04/02/2020 11:30:27.749 AM
%
% Multimodal Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca