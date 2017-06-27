function test_bug1125

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_preproc_bandpassfilter ft_preproc_bandstopfilter ft_preproc_lowpassfilter ft_preproc_highpassfilter

N = 1000;
x1 = randn(1,N)+5*rand(1);

y19 = ft_preproc_bandpassfilter(x1, 1000, [15 25],[],'firls');
y29 = ft_preproc_bandpassfilter(x1, 1000, [15 25],[],'firls');
figure; plot(linspace(0,1000,N),abs(fft(y19))); axis([0 40 0 inf])
hold on; plot(linspace(0,1000,N),abs(fft(y29)),'m'); axis([0 40 0 inf])
% HUGE diff DC

y19 = ft_preproc_bandstopfilter(x1, 1000, [25 15],[],'firls');
y29 = ft_preproc_bandstopfilter(x1, 1000, [25 15],[],'firls');
figure; plot(linspace(0,1000,N),abs(fft(y19))); axis([0 100 0 inf])
hold on; plot(linspace(0,1000,N),abs(fft(y29)),'m'); axis([0 100 0 inf])
figure; plot(y19); axis([0 100 0 inf])
hold on; plot(y29,'m'); axis([0 100 0 inf])
% minor diff DC

y19 = ft_preproc_lowpassfilter(x1, 1000, [15],[],'firls');
y29 = ft_preproc_lowpassfilter(x1, 1000, [15],[],'firls');
figure; plot(linspace(0,1000,N),abs(fft(y19))); axis([0 40 0 inf])
hold on; plot(linspace(0,1000,N),abs(fft(y29)),'m'); axis([0 40 0 inf])
% minor diff DC

y19 = ft_preproc_highpassfilter(x1, 1000, [15],[],'firls');
y29 = ft_preproc_highpassfilter(x1, 1000, [15],[],'firls');
figure; plot(linspace(0,1000,N),abs(fft(y19))); axis([0 40 0 inf])
hold on; plot(linspace(0,1000,N),abs(fft(y29)),'m'); axis([0 40 0 inf])
% HUGE diff DC

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1129.mat

cfg = [];
cfg.bpfreq = [15 25];
cfg.bpfilter = 'yes';
cfg.demean = 'yes';
fdmo = ft_preprocessing(cfg,raw3); %old
fdmn = ft_preprocessing(cfg,raw3); %new

cfg = [];
cfg.bpfreq = [15 25];
cfg.bpfilter = 'yes';
fo = ft_preprocessing(cfg,raw3);
fn = ft_preprocessing(cfg,raw3);

cfg = [];
cfg.lpfreq = [15];
cfg.lpfilter = 'yes';
lo = ft_preprocessing(cfg,raw3);
ln = ft_preprocessing(cfg,raw3);

cfg = [];
cfg.bsfreq = [15 25];
cfg.bsfilter = 'yes';
so = ft_preprocessing(cfg,raw3);
sn = ft_preprocessing(cfg,raw3);

