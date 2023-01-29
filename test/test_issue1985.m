function test_issue1985

% WALLTIME 00:10:00
% MEM 2gb

% DEPENDENCY filter_with_correction fir_filterdcpadded

nchan = 1;
fsample = 1000;
nsample = 2*fsample;

dat = zeros(nchan, nsample);
dat(1,fsample) = 1;

%%

fhp = 0.3;
flp = 30;
fbp = [fhp flp];

order = 30;
type = 'firws';

% use defaulls where possible
filt1 = ft_preproc_bandpassfilter(dat, fsample, fbp, order, type, 'onepass-zerophase');
filt2 = ft_preproc_bandpassfilter(dat, fsample, fbp, order, type, 'onepass-reverse-zerophase');
filt3 = ft_preproc_bandpassfilter(dat, fsample, fbp, order, type, 'onepass-minphase');


%%

figure; hold on
plot(dat * 0.1)
plot(filt1+0.01)
plot(filt2+0.02)
plot(filt3+0.03)
legend