% FT_EXAMPLE_READ_MEF3P0 an example to read MEF 3.0 data set into FieldTrip

%% set the parameters
clear

% look for session path of data
p = mfilename('fullpath');
sesspath = fullfile(fileparts(p), 'mef_3p0.mefd');

% set the password
password = struct('Level1Password', 'password1', 'Level2Password',...
    'password2', 'AccessLevel', 2);

%% read the MEF 3.0 data into MATLAB using class MEFFieldTrip_3p0

% get the object to read data into MATLAB
mef_ft = MEFFieldTrip_3p0(sesspath, password);

% inspect the channels in the session
% note that the channel names are in alphabetic order, which is the default
% order
acq_chan_num = num2str(mef_ft.getAcqChanNumber);
fprintf('Acqusition channel number: %s\n', acq_chan_num)
disp('Channel Names: ')
disp(mef_ft.ChannelName)

% let's sort the channel names according to the number of the acquisition
% channels
mef_ft = MEFFieldTrip_3p0(sesspath, password, 'SortChannel', 'number');
acq_chan_num = num2str(mef_ft.getAcqChanNumber);
fprintf('Acqusition channel number: %s\n', acq_chan_num)
disp('Channel Names: ')
disp(mef_ft.ChannelName)

% now let's import the the first 10 seconds data of two channels 
% 'Left_Occipital-Ref' and 'left-right_occipital' into MATLAB
mef_ft.SelectedChannel = ["Left_Occipital-Ref" "left-right_occipital"];
mef_ft.StartEnd = [0 10]; % you can specify the number of samples too
                          % for examples [1 2561] or seconds [0 10]
mef_ft.SEUnit = 'Second'; % the unit for reading data can be 'Index',
                          % 'Second', 'uUTC', 'Minute', 'Hour' and 'Day'
[X, t] = mef_ft.importSession;
fs = mef_ft.SamplingFrequency;

figure
ph = plot(t/fs, X);
legend(ph, mef_ft.SelectedChannel)
xlim([0, 1]) % zoom into the 1st second data
xlabel('Time (s)')
ylabel('Amplitude')

%% read the MEF 3.0 data using FieldTrip routines
% read data header with ft_read_header()
% --------------------------------------
% sort channel using key-value pair: 'SortChannel' can be either 'alphebat'
% (default) or 'number', which sorts the channel according to the
% acqusition_channel_number.
hdr = ft_read_header(sesspath, 'password', password, 'SortChannel', 'number');

% read a specific channel with ft_read_data()
% Note: if the key-value 'SortChannel' is not provided to ft_read_data(),
% the function uses header.SortChannel. If 'header' is not provided, the
% function use default 'alphabet'.  If both 'header' and 'SortChannel' are
% provided, there will be an warning of the conflict and the function
% chooses header.SortChannel. The same rule is applied to function
% ft_read_event().
chpath = fullfile(sesspath, [hdr.label{4}, '.timd']);
x = ft_read_data(chpath, 'begsample', 1, 'endsample', 2561, 'header', hdr,...
    'password', password, 'chanindx', 1); % don't ommit 'chanindx'
figure
plot(x)
xlim([1 256])
xlabel('Time (sample index)')
legend(hdr.label{4})

%% read data with ft_read_data() but specifying time interval using seconds
% -------------------------------------------------------------------------
% Let's import 10 seconds data at the beginning of the recording with an
% assumption that the trigger was at 0 second.
in_unit = 'second';
be_second = [0, 10, 0]; % 10-second time of data from the start
out_unit = 'index';
be_sample = mef_ft.SessionUnitConvert(be_second, in_unit, out_unit);
dat = ft_read_data(sesspath,...
    'begsample', be_sample(1),...
    'endsample', be_sample(2),...
    'header', hdr,...
    'password', password,...
    'chanindx', [4 1 2 3]); % the order of read data can be decided by
                            % key-value 'chanindx'

t = linspace(be_second(1), be_second(2), be_sample(2)-be_sample(1)+1);
figure
plot(t, dat')
xlim([0 1])
xlabel('Time (s)')
legend(hdr.label{4}, hdr.label{1}, hdr.label{2}, hdr.label{3})


% read data with ft_preprocessing()
% ---------------------------------
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
% Created: Sun 03/22/2020  9:03:27.318 PM
% Revision: 0.5  Date: Fri 04/10/2020 12:33:09.126 AM
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)