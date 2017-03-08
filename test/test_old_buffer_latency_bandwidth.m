function test_old_buffer_latency_bandwidth

% MEM 1gb
% WALLTIME 00:10:00


is_octave=~ft_platform_supports('matlabversion',1,inf);
if is_octave
  % When using Octave, /realtime/src/buffer/matlab/buffer.c gave
  % problems compiling.
  % The error raised when compiling is:
  %   In file included from buffer.c:11:
  %   In file included from $OCTAVE/4.0.3/include/octave-4.0.3/octave/matrix.h:30:
  %   In file included from $OCTAVE/4.0.3/include/octave-4.0.3/octave/mx-base.h:28:
  %   $OCTAVE/4.0.3/include/octave-4.0.3/octave/MatrixType.h:27:1: error:
  %     unknown type name 'class'
  %   class Matrix;
  %   ^
  %
  % skip the test
  reason=sprintf(['%s requires a compiled version of the file '...
                '/realtime/src/buffer/matlab/buffer.c, but in Octave '...
                'attempts to compile it have not been succesful. '...
                'Therefore this test has been disabled'],...
                mfilename());
  moxunit_throw_test_skipped_exception(reason);
end

filename = 'buffer://localhost:1972';

hdr = [];
hdr.nChans = 10;
hdr.nSamples = 0;
hdr.nSamplesPre = 0;
hdr.Fs = 1000;
% hdr.label = ...

ft_write_data(filename, [], 'header', hdr, 'append', 0);

siz = round(logspace(0, 3, 100));
dat = zeros(hdr.nChans, max(siz));
stopwatch = tic;
t(1) = toc(stopwatch);
for i=1:length(siz)
  disp(i)
  ft_write_data(filename, dat(:,1:siz(i)), 'header', hdr, 'append', 1);
  t(i+1) = toc(stopwatch);
end

t = diff(t);
siz = siz*hdr.nChans;

figure
plot(siz, t, '.');

[p, s] = polyfit(siz, t, 1);
sel = round(0.33*numel(siz)):numel(siz);
[p, s] = polyfit(siz(sel), t(sel), 1);
e = polyval(p, siz);

hold on
plot(siz, e, 'r-')

fprintf('latency   = %f sec\n', p(2));
fprintf('bandwidth = %f samples/sec\n', 1/p(1));

