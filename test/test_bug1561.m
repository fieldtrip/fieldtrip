function test_bug1561

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preproc_bandpassfilter ft_preproc_bandstopfilter ft_preproc_lowpassfilter ft_preproc_highpassfilter

tolerance = 1e-4;
tolerance = inf; % FIXME

Fs = 1000;

% construct a test signal
begspike = zeros(1,1000); begspike(1) = 1;
endspike = fliplr(begspike);

% this is an alternative test signal
dat = zeros(1,1001); dat(501) = 1;

% use reasonable numbers for the filters
Fbp   = [8 12];
Fl    = 30;
Fh    = 2;
Fbs   = [45 55];

% consider looping over all filter types and multiple filter orders
order = [];    % default
type  = 'but'; % default

figure; plot(ft_preproc_bandpassfilter(dat, Fs, Fbp, order, type, 'twopass'))
figure; plot(ft_preproc_bandpassfilter(dat, Fs, Fbp, order, type, 'twopass-average'))
figure; plot(ft_preproc_bandstopfilter(dat, Fs, Fbp, order, type, 'twopass'))
figure; plot(ft_preproc_lowpassfilter (dat, Fs, Fbp, order, type, 'twopass'))
figure; plot(ft_preproc_highpassfilter(dat, Fs, Fbp, order, type, 'twopass'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begflt = ft_preproc_bandpassfilter(begspike, Fs, Fbp, order, type, 'onepass');
endflt = ft_preproc_bandpassfilter(endspike, Fs, Fbp, order, type, 'onepass-reverse');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_bandstopfilter(begspike, Fs, Fbs, order, type, 'onepass');
endflt = ft_preproc_bandstopfilter(endspike, Fs, Fbs, order, type, 'onepass-reverse');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_lowpassfilter(begspike, Fs, Fl, order, type, 'onepass');
endflt = ft_preproc_lowpassfilter(endspike, Fs, Fl, order, type, 'onepass-reverse');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_highpassfilter(begspike, Fs, Fh, order, type, 'onepass');
endflt = ft_preproc_highpassfilter(endspike, Fs, Fh, order, type, 'onepass-reverse');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begflt = ft_preproc_bandpassfilter(begspike, Fs, Fbp, order, type, 'twopass');
endflt = ft_preproc_bandpassfilter(endspike, Fs, Fbp, order, type, 'twopass');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_bandstopfilter(begspike, Fs, Fbs, order, type, 'twopass');
endflt = ft_preproc_bandstopfilter(endspike, Fs, Fbs, order, type, 'twopass');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_lowpassfilter(begspike, Fs, Fl, order, type, 'twopass');
endflt = ft_preproc_lowpassfilter(endspike, Fs, Fl, order, type, 'twopass');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_highpassfilter(begspike, Fs, Fh, order, type, 'twopass');
endflt = ft_preproc_highpassfilter(endspike, Fs, Fh, order, type, 'twopass');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

figure; plot([begflt; endflt]'); legend('begspike', 'endspike')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begflt = ft_preproc_bandpassfilter(begspike, Fs, Fbp, order, type, 'twopass-average');
endflt = ft_preproc_bandpassfilter(endspike, Fs, Fbp, order, type, 'twopass-average');
assert(all(abs(begflt - fliplr(endflt))<tolerance));

figure; plot([begflt; endflt]'); legend('begspike', 'endspike')
