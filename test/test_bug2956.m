function test_bug2956

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_preprocessing preproc ft_preproc_dftfilter 
% Test Script for the new dftfilter option "Spectrum Interpolation" to reduce power line noise in data

%% Simulate some data with power line noise: 50 Hz sinusoid (increasing amplitude) added to noise
% (or load data that contains line noise)

len_sec = 100;              % length of signal in seconds
Fs = 1000;                  % SR
Tm = 1/Fs;  
Ln = Fs.*len_sec;

amp = (0.1.*(1:Ln))/Fs;      % increasing amplitude values 
timepts = (0:Ln-1)*Tm;               
PLN_50Hz = amp.*sin(2*pi*50*timepts);    
PLN_100Hz = amp/2.*sin(2*pi*100*timepts); 
PLN_150Hz = amp/3.*sin(2*pi*150*timepts);
s1 = PLN_50Hz + PLN_100Hz + PLN_150Hz + randn(size(timepts)); % random noise + 50 Hz power line noise (+harmonics)  

% add S1 to dummy data structure
data = [];
data.trial{1} = s1;
data.fsample = Fs;
data.label = {'RandomNoise_&_PowerLineNoise'};
data.time{1} = timepts;
data.fsample = Fs;
data.sampleinfo = rand(1,2);

% plot simulated data
figure; plot(timepts(1:Ln),data.trial{1}(1:Ln));

% plot power spectrum of simulated data 
data_fft = fft(data.trial{1},Ln,2); 
frq = Fs*linspace(0,1,Ln+1);
figure; 
semilogy(frq(1:Ln),abs(data_fft(1,1:Ln).^2));
freq_res = frq(2); 
    
%% apply "Spectrum Interpolation" to reduce 50 Hz noise as new dft filter method 

cfg= [];
cfg.dftfilter = 'yes';
cfg.dftfreq = [50 100 150]; 
cfg.dftreplace = 'neighbour'; % implicates spectrum interpolation
cfg.dftbandwidth = [2 2 2]; 
cfg.dftneighbourwidth = [2 2 2]; 
data_intpl = ft_preprocessing(cfg, data);

% plot simulated data after spectrum interpolation, superimposed on original noisy data
figure; plot(timepts, data.trial{1},'red');
hold; plot(timepts,data_intpl.trial{1});


%  plot power spectrum after spectrum interpolation
data_fft2 = fft(data_intpl.trial{1},Ln,2);
figure;semilogy(frq(1:Ln),abs(data_fft2(1,1:Ln)).^2);

%% now try "DFT filter" option, doesn't work for varying 50 Hz amplitude, as shown in the fieldtrip tutorial

cfg= [];
cfg.dftfilter = 'yes';
cfg.dftfreq = [50 100 150]; %
cfg.dftreplace = 'zero'; % implies dft filter (formerly the only option)
data_dftfilt = ft_preprocessing(cfg, data);

% plot simulated data after dft filtering
figure;  plot(timepts, data.trial{1},'red');
hold; plot(timepts,data_dftfilt.trial{1});

%  plot power spectrum after dft filtering
data_fft2 = fft(data_dftfilt.trial{1},Ln,2);
figure;semilogy(frq(1:Ln),abs(data_fft2(1,1:Ln)).^2);

%%
