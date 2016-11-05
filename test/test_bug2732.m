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

warning('');          % clear previous warnings
ft_warning('-clear'); % clear previous warnings

%%

out0 = ft_preproc_baselinecorrect(dat0);
out1 = ft_preproc_baselinecorrect(dat1);
out2 = ft_preproc_baselinecorrect(dat2);
out3 = ft_preproc_baselinecorrect(dat3);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

Fbp = [5 15];
out0 = ft_preproc_bandpassfilter(dat0,Fs,Fbp);
out1 = ft_preproc_bandpassfilter(dat1,Fs,Fbp);
out2 = ft_preproc_bandpassfilter(dat2,Fs,Fbp);
out3 = ft_preproc_bandpassfilter(dat3,Fs,Fbp);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

Fbp = [5 15];
out0 = ft_preproc_bandstopfilter(dat0, Fs, Fbp);
out1 = ft_preproc_bandstopfilter(dat1, Fs, Fbp);
out2 = ft_preproc_bandstopfilter(dat2, Fs, Fbp);
out3 = ft_preproc_bandstopfilter(dat3, Fs, Fbp);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

Fhp = 5;
out0 = ft_preproc_highpassfilter(dat0, Fs, Fhp);
out1 = ft_preproc_highpassfilter(dat1, Fs, Fhp);
out2 = ft_preproc_highpassfilter(dat2, Fs, Fhp);
out3 = ft_preproc_highpassfilter(dat3, Fs, Fhp);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

Flp = 15;
out0 = ft_preproc_lowpassfilter(dat0, Fs, Flp);
out1 = ft_preproc_lowpassfilter(dat1, Fs, Flp);
out2 = ft_preproc_lowpassfilter(dat2, Fs, Flp);
out3 = ft_preproc_lowpassfilter(dat3, Fs, Flp);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

order = 10;
out0 = ft_preproc_medianfilter(dat0, order);
out1 = ft_preproc_medianfilter(dat1, order);
out2 = ft_preproc_medianfilter(dat2, order);
out3 = ft_preproc_medianfilter(dat3, order);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));

assert(sum(isnan(out1(:,end)))==0);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==order);

%%

Fl = 50;
out0 = ft_preproc_dftfilter(dat0, Fs, Fl);
out1 = ft_preproc_dftfilter(dat1, Fs, Fl);
out2 = ft_preproc_dftfilter(dat2, Fs, Fl);
out3 = ft_preproc_dftfilter(dat3, Fs, Fl);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

refdat = randn(2,ntime);
out0 = ft_preproc_denoise(dat0, refdat);
out1 = ft_preproc_denoise(dat1, refdat);
out2 = ft_preproc_denoise(dat2, refdat);
out3 = ft_preproc_denoise(dat3, refdat);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

out0 = ft_preproc_derivative(dat0);
out1 = ft_preproc_derivative(dat1);
out2 = ft_preproc_derivative(dat2);
out3 = ft_preproc_derivative(dat3);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==0);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==0);
assert(sum(isnan(out3(10,:)))==2);

%%

out0 = ft_preproc_detrend(dat0);
out1 = ft_preproc_detrend(dat1);
out2 = ft_preproc_detrend(dat2);
out3 = ft_preproc_detrend(dat3);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

out0 = ft_preproc_hilbert(dat0);
out1 = ft_preproc_hilbert(dat1);
out2 = ft_preproc_hilbert(dat2);
out3 = ft_preproc_hilbert(dat3);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

order = 2;
out0 = ft_preproc_polyremoval(dat0, order);
out1 = ft_preproc_polyremoval(dat1, order);
out2 = ft_preproc_polyremoval(dat2, order);
out3 = ft_preproc_polyremoval(dat3, order);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);

%%

out0 = ft_preproc_rereference(dat0);
out1 = ft_preproc_rereference(dat1);
out2 = ft_preproc_rereference(dat2);
out3 = ft_preproc_rereference(dat3);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==0);
assert(sum(isnan(out2(:,end)))==nchan);
assert(sum(isnan(out3(:,10)))==nchan);
assert(sum(isnan(out3(10,:)))==1);

%%

Fnew = 500;
method = 'resample';
out0 = ft_preproc_resample(dat0, Fs, Fnew, method);
out1 = ft_preproc_resample(dat1, Fs, Fnew, method);
out2 = ft_preproc_resample(dat2, Fs, Fnew, method);
out3 = ft_preproc_resample(dat3, Fs, Fnew, method);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==0);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==15); % the NaN only has a local effect

%%

width = 5;
out0 = ft_preproc_slidingrange(dat0, width);
out1 = ft_preproc_slidingrange(dat1, width);
out2 = ft_preproc_slidingrange(dat2, width);
out3 = ft_preproc_slidingrange(dat3, width);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==0);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==0);
assert(sum(isnan(out3(10,:)))==0);

%%

width = 5;
out0 = ft_preproc_smooth(dat0, width);
out1 = ft_preproc_smooth(dat1, width);
out2 = ft_preproc_smooth(dat2, width);
out3 = ft_preproc_smooth(dat3, width);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==0);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==width);

%%

out0 = ft_preproc_standardize(dat0);
out1 = ft_preproc_standardize(dat1);
out2 = ft_preproc_standardize(dat2);
out3 = ft_preproc_standardize(dat3);

[msg, id] = lastwarn;
assert(isequal(id, 'FieldTrip:dataContainsNaN'));
warning(''); % clear previous warnings

assert(sum(isnan(out1(:,end)))==nchan);
assert(sum(isnan(out2(:,end)))==1);
assert(sum(isnan(out3(:,10)))==1);
assert(sum(isnan(out3(10,:)))==ntime);


