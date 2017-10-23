function test_bug1770

% MEM 1500mb
% WALLTIME 00:10:00

% TEST read_neuralynx_dma read_neuralynx_ncs

filenameA = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1770/2012-07-14_15-33-09');
filenameB = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1770/DigitalLynxRawDataFile.nrd');

hdrA = ft_read_header(filenameA);
hdrB = ft_read_header(filenameB);

begsample = 1;
endsample = 100000;

datA = ft_read_data(filenameA, 'begsample', begsample, 'endsample', endsample);
datB = ft_read_data(filenameB, 'begsample', begsample, 'endsample', endsample);

ch1A = double(datA( 1,:));
ch1B = double(datB(18,:));

figure
subplot(2,1,1);
plot(ft_preproc_standardize(ch1A))
subplot(2,1,2);
plot(-ft_preproc_standardize(ch1B))

ch1As =  ft_preproc_standardize(ch1A);
ch1Bs = -ft_preproc_standardize(ch1B);

figure
hold on
plot(ch1As, 'b.');
plot(ch1Bs, 'r.');

% sofar there seems no reason for concern w.r.t. the reading function
% let us look at the filter kernel and the delay

hA = fft(ch1As);
hB = fft(ch1Bs);
H = hA./hB;
figure
plot(ifft(H));

