function test_old_preproc_resample

% MEM 1gb
% WALLTIME 00:10:00


%this script tests the anti-aliasing filter in the resample function as used by ft_preproc_resample

fs    = 1000;
nsec  = 100;
nsmp  = fs*nsec;
tim   = ([1:nsmp]-1)./fs;
foi   = [0:(fs-1)];
 
datin       = randn(1,nsmp);
datin2      = reshape(datin, [fs nsec]);
datin2pow   = abs(fft(datin2, [], 1)).^2;

dathigh     = ft_preproc_highpassfilter(datin, fs, 300, 10, 'but', 'twopass');
dathigh     = dathigh + 0.1.*sin(2.*pi.* 400 .* tim);
dathigh2    = reshape(dathigh, [fs nsec]);
dathigh2pow = abs(fft(dathigh2, [], 1)).^2;

datlow     = ft_preproc_lowpassfilter(datin, fs, 300, 10, 'but', 'twopass');
datlow     = datlow + 0.1.*sin(2.*pi.* 400 .* tim);
datlow2    = reshape(datlow, [fs nsec]);
datlow2pow = abs(fft(datlow2, [], 1)).^2;

figure;hold on;
plot(foi, mean(datin2pow,2));
plot(foi, mean(dathigh2pow,2),'r');
plot(foi, mean(datlow2pow,2), 'm');

dathigh     = ft_preproc_highpassfilter(datin, fs, 400, 10, 'fir', 'twopass');
dathigh     = dathigh + 0.1.*sin(2.*pi.* 350 .* tim) + 0.1.*randn(1,nsmp);
dathigh2    = reshape(dathigh, [fs nsec]);
dathigh2pow = abs(fft(dathigh2, [], 1)).^2;

datlow     = ft_preproc_lowpassfilter(datin, fs, 400, 10, 'fir', 'twopass');
datlow     = datlow + 0.1.*sin(2.*pi.* 350 .* tim);
datlow2    = reshape(datlow, [fs nsec]);
datlow2pow = abs(fft(datlow2, [], 1)).^2;

figure;hold on;
plot(foi, mean(datin2pow,2));
plot(foi, mean(dathigh2pow,2),'r');
plot(foi, mean(datlow2pow,2), 'm');

newfs    = 500;
newfoi   = [0:(newfs-1)];
dathighr = ft_preproc_resample(dathigh, fs, newfs, 'resample');
dathighx = reshape(dathighr, [newfs nsec]);
dathighxpow = abs(fft(dathighx, [], 1)).^2;

dathighlow = ft_preproc_resample(ft_preproc_lowpassfilter(dathigh, fs, 150, 2, 'but', 'twopass'), fs, newfs, 'resample');
dathighlow = reshape(dathighlow, [newfs nsec]);
dathighlowpow = abs(fft(dathighlow, [], 1)).^2;

figure;hold on;
plot(foi(1:250),    (mean(dathigh2pow(1:250,2:end-1),2)));
plot(newfoi(1:250), (mean(dathighxpow(1:250,2:end-1),2)),'r');
plot(newfoi(1:250), (mean(dathighlowpow(1:250,2:end-1),2)),'m');

figure;hold on;
plot(foi(1:250), mean(dathighxpow(1:250,2:end-1),2)./mean(dathigh2pow(1:250,2:end-1),2));
