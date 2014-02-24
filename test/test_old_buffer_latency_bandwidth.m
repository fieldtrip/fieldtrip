function test_old_buffer_latency_bandwidth

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_buffer_latency_bandwidth

filename = 'buffer://localhost:1972';

hdr = [];
hdr.nChans = 10;
hdr.nSamples = 0;
hdr.nSamplesPre = 0;
hdr.Fs = 1000;
% hdr.label = ...

write_data(filename, [], 'header', hdr, 'append', 0);

siz = round(logspace(0, 3, 100));
dat = zeros(hdr.nChans, max(siz));
stopwatch = tic;
t(1) = toc(stopwatch);
for i=1:length(siz)
disp(i)
write_data(filename, dat(:,1:siz(i)), 'header', hdr, 'append', 1);
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

