function test_bug2732

% WALLTIME 00:10:00
% MEM 1gb

nchan = 16;
ntime = 1000;
Fs = 1000;

dat0 = randn(nchan, ntime);

dat1 = dat0; dat1( :,10) = nan;
dat2 = dat0; dat2(10, :) = nan;
dat3 = dat0; dat3(10,10) = nan;

%%

out0 = ft_preproc_baselinecorrect(dat0);
out1 = ft_preproc_baselinecorrect(dat1);
out2 = ft_preproc_baselinecorrect(dat2);
out3 = ft_preproc_baselinecorrect(dat3);

%%

Fbp = [5 15];
out0 = ft_preproc_bandpassfilter(dat0,Fs,Fbp);
out1 = ft_preproc_bandpassfilter(dat1,Fs,Fbp);
out2 = ft_preproc_bandpassfilter(dat2,Fs,Fbp);
out3 = ft_preproc_bandpassfilter(dat3,Fs,Fbp);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));

%%

Fbp = [5 15];
out0 = ft_preproc_bandstopfilter(dat0,Fs,Fbp);
out1 = ft_preproc_bandstopfilter(dat1,Fs,Fbp);
out2 = ft_preproc_bandstopfilter(dat2,Fs,Fbp);
out3 = ft_preproc_bandstopfilter(dat3,Fs,Fbp);

%%

ft_preproc_highpassfilter(dat);
ft_preproc_lowpassfilter(dat);
ft_preproc_medianfilter(dat);

ft_preproc_denoise(dat);
ft_preproc_derivative(dat);
ft_preproc_detrend(dat);
ft_preproc_dftfilter(dat);
ft_preproc_hilbert(dat);
ft_preproc_padding(dat);
ft_preproc_polyremoval(dat);
ft_preproc_rectify(dat);
ft_preproc_rereference(dat);
ft_preproc_resample(dat);
ft_preproc_slidingrange(dat);
ft_preproc_smooth(dat);
ft_preproc_standardize(dat);

