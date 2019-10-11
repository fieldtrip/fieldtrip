% $Id: realign_synth_pca.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

%% build synthetic data

% synthetic data
randn('seed',0);
rand('seed',0);

clear EEG

EEG.pnts = 500;
EEG.srate = 1000; % 0.5 sec. recording
EEG.trials = 100;
EEG.nbchan = 1;

EEG.times = linspace(0,EEG.pnts / EEG.srate,EEG.pnts); % tf*100

noiselevel = 5;
[EEG,ref,latencies,snr] = generate_synth_eeg_data(EEG,noiselevel);

X = squeeze(EEG.data)';
X = X - ones(size(X,1),1)*mean(X,1);
norms = sqrt(sum(X.^2,2));

X = diag(1 ./ norms) * X;
[PC, SCORE, LATENT]=princomp(X);

[tmp,order] = sort(SCORE(:,1));

data = squeeze(EEG.data)';

%% view results
figure(400)
plot(SCORE(:,1),SCORE(:,2),'+');
xlabel('X');
ylabel('Y');
% title('Points in 2D with PCA');
% savefig('pca_embedding')

figure(401)
imagesc(data,[min(data(:)),max(data(:))])
colorbar
xlabel('Time (ms)')
ylabel('Trial')
title('Data')
set(gca,'XTick',[])

figure(402)
imagesc(data(order,:),[min(data(:)),max(data(:))])
colorbar
xlabel('Time (ms)')
ylabel('Trial')
% title('Data Ordered with PCA')
% set(gcf,'WindowButtonDownFcn','callbackRaster(EEG)'); % set the callback
set(gca,'XTick',[])
% savefig('pca_realigned')


