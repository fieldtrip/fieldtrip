function test_external_signal

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY hanning

filelist = {'barthannwin'
            'blackmanharris'
            'bohmanwin'
            'boxcar'
            'butter'
            'filtfilt'
            'flattopwin'
            'gausswin'
            'hann'
            'hanning'
            'hilbert'
            'kaiser'
            'nuttallwin'
            'parzenwin'
            'rectwin'
            'tukeywin'
            'triang'
            'window'};

[ftver, ftpath] = ft_version;
restoredefaultpath

% ensure that the the 'compat' (i.e. external/signal) is used
global ft_default
ft_default.toolbox.signal = 'compat';
addpath(ftpath);
ft_defaults;

for k = 1:numel(filelist)
  assert(exist(filelist{k}, 'file'))
  switch filelist{k}
    case 'barthannwin'
      assert(isequal(barthannwin(1), 1));
      assert(isequal(barthannwin (2), zeros(2, 1)));
    case 'blackmanharris'
      assert(isequal(blackmanharris(1),  1));
      assert(isalmostequal(blackmanharris(2),  0.00006 * ones(2, 1), 'abstol', eps));
      assert(isalmostequal(blackmanharris(15), flipud(blackmanharris(15)), 'abstol', 10*eps));
      assert(isalmostequal(blackmanharris(16), flipud(blackmanharris(16)), 'abstol', 10*eps));
      assert(isequal(blackmanharris(15), blackmanharris(15, 'symmetric')));
      tmp = blackmanharris(16);
      assert(isequal(tmp(1:15), blackmanharris(15, 'periodic')));
    case 'bohmanwin'
%assert(isequal(bohmanwin (1), 1)
%assert(isequal(bohmanwin (2), zeros (2, 1))
    case 'boxcar'
%assert(isequal(boxcar (1), 1)
%assert(isequal(boxcar (2), ones (2, 1))
%assert(isequal(boxcar (100), ones (100, 1))
    case 'butter'
%!shared sf, sf2, off_db
%! off_db = 0.5;
%! ## Sampling frequency must be that high to make the low pass filters pass.
%! sf = 6000; sf2 = sf/2;
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];

%!test
%! ## Test low pass order 1 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 1, 50 / sf2 );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 4 ) - damp_db( 5 ), damp_db( 1 : 3 ) ], [ 6 0 0 -3 ], off_db )

%!test
%! ## Test low pass order 4 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 4, 50 / sf2 );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 4 ) - damp_db( 5 ), damp_db( 1 : 3 ) ], [ 24 0 0 -3 ], off_db )

%!test
%! ## Test high pass order 1 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 1, 50 / sf2, "high" );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 2 ) - damp_db( 1 ), damp_db( 3 : end ) ], [ 6 -3 0 0 ], off_db )

%!test
%! ## Test high pass order 4 with 3dB @ 50Hz
%! data=[sinetone(5,sf,10,1),sinetone(10,sf,10,1),sinetone(50,sf,10,1),sinetone(200,sf,10,1),sinetone(400,sf,10,1)];
%! [b, a] = butter ( 4, 50 / sf2, "high" );
%! filtered = filter ( b, a, data );
%! damp_db = 20 * log10 ( max ( filtered ( end - sf : end, : ) ) );
%! assert ( [ damp_db( 2 ) - damp_db( 1 ), damp_db( 3 : end ) ], [ 24 -3 0 0 ], off_db )

%% Test input validation
%!error [a, b] = butter ()
%!error [a, b] = butter (1)
%!error [a, b] = butter (1, 2, 3, 4, 5)
%!error [a, b] = butter (.5, .2)
%!error [a, b] = butter (3, .2, "invalid")

%!error [a, b] = butter (9, .6, "stop")
%!error [a, b] = butter (9, .6, "bandpass")

%!error [a, b] = butter (9, .6, "s", "high")

%% Test output orientation
%!test
%! butter (9, .6);
%! assert (isrow (ans));
%!test
%! A = butter (9, .6);
%! assert (isrow (A));
%!test
%! [A, B] = butter (9, .6);
%! assert (isrow (A));
%! assert (isrow (B));
%!test
%! [z, p, g] = butter (9, .6);
%! assert (iscolumn (z));
%! assert (iscolumn (p));
%! assert (isscalar (g));
%!test
%! [a, b, c, d] = butter (9, .6);
%! assert (ismatrix (a));
%! assert (iscolumn (b));
%! assert (isrow (c));
%! assert (isscalar (d));
    case 'filtfilt'
%!test
%! randn('state',0);
%! r = randn(1,200);
%! [b,a] = butter(10, [.2, .25]);
%! yfb = filtfilt(b, a, r);
%! assert (size(r), size(yfb));
%! assert (mean(abs(yfb)) < 1e3);
%! assert (mean(abs(yfb)) < mean(abs(r)));
%! ybf = fliplr(filtfilt(b, a, fliplr(r)));
%! assert (mean(abs(ybf)) < 1e3);
%! assert (mean(abs(ybf)) < mean(abs(r)));

%!test
%! randn('state',0);
%! r = randn(1,1000);
%! s = 10 * sin(pi * 4e-2 * (1:length(r)));
%! [b,a] = cheby1(2, .5, [4e-4 8e-2]);
%! y = filtfilt(b, a, r+s);
%! assert (size(r), size(y));
%! assert (mean(abs(y)) < 1e3);
%! assert (corr(s(250:750), y(250:750)) > .95)
%! [b,a] = butter(2, [4e-4 8e-2]);
%! yb = filtfilt(b, a, r+s);
%! assert (mean(abs(yb)) < 1e3);
%! assert (corr(y, yb) > .99)

%!test
%! randn('state',0);
%! r = randn(1,1000);
%! s = 10 * sin(pi * 4e-2 * (1:length(r)));
%! [b,a] = butter(2, [4e-4 8e-2]);
%! y = filtfilt(b, a, [r.' s.']);
%! yr = filtfilt(b, a, r);
%! ys = filtfilt(b, a, s);
%! assert (y, [yr.' ys.']);
%! y2 = filtfilt(b.', a.', [r.' s.']);
%! assert (y, y2);
    case 'flattopwin'
      assert(isequal(flattopwin(1), 1));
      assert(isalmostequal(flattopwin(2), 0.0042 / 4.6402 * ones (2, 1), 'abstol', eps));
      assert(isalmostequal(flattopwin(15), flipud(flattopwin(15)), 'abstol', 10*eps));
      assert(isalmostequal(flattopwin(16), flipud(flattopwin(16)), 'abstol', 10*eps));
      assert(isequal(flattopwin(15), flattopwin(15, 'symmetric')));
      tmp = flattopwin(16);
      assert(isequal(tmp(1:15), flattopwin(15, 'periodic')));
    case 'gausswin'
      assert(isequal(gausswin(1), 1));
      assert(isequal(gausswin(2), [exp(-3.125); exp(-3.125)]));
      assert(isequal(gausswin(3), [exp(-3.125); 1; exp(-3.125)]));
    case 'hann'
      assert(isequal(hann(1), 1));
      assert(isalmostequal(hann(2), [0.75;0.75], 'abstol', eps)); %the original octave required [0;0], yet matlab returns [0.75;0.75]
      assert(isalmostequal(hann(16), flipud(hann(16)), 'abstol', 10*eps));
      assert(isalmostequal(hann(15), flipud(hann(15)), 'abstol', 10*eps));
      assert(isequal(hann(15), hann(15, 'symmetric')));
      tmp = hann(16, 'periodic');
      assert(isequal(tmp(2:16), hann(15, 'symmetric'))); % this is different from the octave test, since it is how matlab does it
    case 'hanning'
    case 'hilbert'
      % FIXME
    case 'kaiser'
      assert(isequal(kaiser(1), 1));
    case 'nuttallwin'
      assert(isequal(nuttallwin(1), 1));
      assert(isalmostequal(nuttallwin(2), zeros(2, 1), 'abstol', eps));
      assert(isalmostequal(nuttallwin(15), flipud(nuttallwin(15)), 'abstol', 10*eps));
      assert(isalmostequal(nuttallwin(16), flipud(nuttallwin(16)), 'abstol', 10*eps));
      assert(isequal(nuttallwin(15), nuttallwin(15, 'symmetric')));
      tmp = nuttallwin(16);
      assert(isequal(tmp(1:15), nuttallwin(15, 'periodic')));
    case 'parzenwin'
      assert(isequal(parzenwin(1), 1));
      assert(isequal(parzenwin(2), 0.25 * ones (2, 1)));
    case 'rectwin'
      assert(isequal(rectwin(1),   1));
      assert(isequal(rectwin(2),   ones(2, 1)));
      assert(isequal(rectwin(100), ones(100, 1)));
    case 'triang'
      assert(isequal(triang(1), 1));
      assert(isequal(triang(2), [1; 1]/2));
      assert(isequal(triang(3), [1; 2; 1]/2));
      assert(isequal(triang(4), [1; 3; 3; 1]/4));
    case 'tukeywin'
      assert(isequal(tukeywin(1), 1));
      assert(isequal(tukeywin(2), zeros (2, 1)));
      assert(isequal(tukeywin(3), [0; 1; 0]));
      assert(isequal(tukeywin(16, 0), rectwin (16)));
      assert(isequal(tukeywin(16, 1), hanning (16)));
    case 'window'
      assert(isequal(window(@hanning, 16), window('hanning', 16)));
      assert(isequal(window(@triang, 16),  window('triang', 16)));
    otherwise
      ft_error('function %s is not part of the official external/signal directory', filelist{k});
  end
end
