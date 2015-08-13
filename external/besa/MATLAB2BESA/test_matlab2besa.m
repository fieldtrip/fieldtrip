%% Examples for saving data matrix as ASCII-vectorized file.

% Load example data 1
load data_avr.mat

% Set parameters
custom_path = './'; % current directory
file_name = 'test1.avr';
data_matrix = data_avr.Data;
time_samples = data_avr.Time;
channel_labels = data_avr.ChannelLabels;
data_scale_factor = 1.0;
time_scale_factor = 1.0;

% Save the data
besa_save2Avr(custom_path, file_name, data_matrix, ...
    time_samples, channel_labels, data_scale_factor, time_scale_factor);

%% Save as generic binary

% Load example data 1
load data_avr.mat

% Set parameters
file_name = 'test1';
data_matrix = data_avr.Data;
SamplingRate = 1000.0/data_avr.DI;

% Save the data
besa_matrix2Gen(data_matrix, SamplingRate, file_name);

%% Load example data 2
load besa_channels.mat

% Set parameters
custom_path = './'; % current directory
file_name = 'test1.avr';
data_matrix = besa_channels.data.amplitudes;
time_samples = besa_channels.data.latencies;
channel_labels = besa_channels.channellabels;
data_scale_factor = 1.0;
time_scale_factor = 1.0;

% Save the data
besa_save2Avr(custom_path, file_name, data_matrix, ...
    time_samples, channel_labels, data_scale_factor, time_scale_factor);