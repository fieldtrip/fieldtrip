% function test_bug1561

% TEST test_bug1561
% TEST ft_preproc_bandpassfilter ft_preproc_lowpassfilter ft_preproc_highpassfilter ft_preproc_bandstopfilter

tolerance = 100*eps;

Fs = 1000;
begspike = zeros(1,1000); begspike(1)   = 1;
endspike = fliplr(begspike);

% use reasonable numbers for the filters
Fbp   = [8 12];
Fl    = 30;
Fh    = 2;
Fbs   = [45 55];

% consider looping over all filter types and multiple filter orders
order = [];    % default
type  = 'but'; % default

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

plot([begflt; endflt]'); legend('begspike', 'endspike')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

begflt = ft_preproc_bandpassfilter(begspike, Fs, Fbp, order, type, 'twopass');
endflt = ft_preproc_bandpassfilter(endspike, Fs, Fbp, order, type, 'twopass');
% assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_bandstopfilter(begspike, Fs, Fbs, order, type, 'twopass');
endflt = ft_preproc_bandstopfilter(endspike, Fs, Fbs, order, type, 'twopass');
% assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_lowpassfilter(begspike, Fs, Fl, order, type, 'twopass');
endflt = ft_preproc_lowpassfilter(endspike, Fs, Fl, order, type, 'twopass');
% assert(all(abs(begflt - fliplr(endflt))<tolerance));

begflt = ft_preproc_highpassfilter(begspike, Fs, Fh, order, type, 'twopass');
endflt = ft_preproc_highpassfilter(endspike, Fs, Fh, order, type, 'twopass');
% assert(all(abs(begflt - fliplr(endflt))<tolerance));

