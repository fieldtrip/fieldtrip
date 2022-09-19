% This script provides a tutorial application of the Gaussian-Copula Mutual
% Information (GCMI) estimator  with a continuous auditory stimulus feature 
% (low pass filtered  speech envelope) in a continuous design

% This script uses cell mode - cells are delimited by %% lines and can 
% be run with: 
% ctrl-enter (windows, linux) or cmd-enter (mac)
% or "Run Section" button from the toolbar 
% or right click -> Evaluate currect section

% You need the gcmi/matlab directory on your matlab path.
% If you are running from the gcmi project tree you can add this with the
% following command:
% addpath('../matlab')

% Questions / comments : robince@gmail.com

%% Download and load data

% this will attempt to download the eeg data from the internet (~86MB)
% alternatively you can download the meg_speech.mat manually
% and place it in the same directory as this script
fname = 'meg_speech.mat';
data_url = ['https://www.robince.net/data/gcmi/data/' fname];

% checks for data in current working directory
if ~exist(fname)
    disp(sprintf('Downloading %s ...', data_url))
    websave(fname,data_url);
    disp('Done.')
end

load(fname);
% csddat : [channels x time]
[Nch, Nt] = size(plndat);
Fs = 50;
% 2 directional components for each channel
Nch = Nch/2;
% permute to have samples first axis
plnH = permute(plndat(1:Nch,:), [2 1]);
plnV = permute(plndat((Nch+1):end,:), [2 1]);

% sum of square of compoments (vector amplitude)
plnsum = sqrt(plnH.^2 + plnV.^2);

% fltspc : low-pass filtered speech envelope of auditory stimulus
% delays to consider in samples (50Hz)
delay_samples = 1:17;
delays = delay_samples * 20; % ms
Ndel = length(delay_samples);

%% Calculate GCMI at a specific sensor and delay

% this is the strongest effect
di = 5;
chi = 231;

% we delay the MEG signal with respect to the stimulus
d = delay_samples(di);
dspeech = fltspc(1:(end-d));
dmeg = plnsum((1+d):end, chi);

% we can then correlate these signals
p = corr(dmeg, dspeech, 'type', 'Pearson')
r = corr(dmeg, dspeech, 'type', 'Spearman')

% we call the GCMI function in a similar fashion to corr
I = gcmi_cc(dmeg, dspeech)

%% Calculate GCMI across all sensors and delays

% following the commonly used mass-univeriate approach we repeat 
% the calculation above for each sensor and delay
% we could call gcmi_cc as above inside the loop, but here we first
% normalise the data separately so we can reuse it for the permutation
% testing

% Gaussian-copula normalisation, works along first axis, applied to each 
% other dimension independently
cplnsum = copnorm(plnsum);
cspeech = copnorm(fltspc);

Iplnsum = zeros(Nch,Ndel);
for di=1:Ndel
    d = delay_samples(di);
    % resize speech according to the lag/delay considered
    dspeech = cspeech(1:(end-d));
    for chi=1:Nch
        % as the data has been copula-normalised we can use the 
        % Gaussian parametric estimator (whatever the distribution was
        % originally)
        Iplnsum(chi,di) = mi_gg(cplnsum((1+d):end,chi), dspeech, true, true);
    end
end

figure
imagesc(delays,[],Iplnsum)
xlabel('MEG-Stimulus Delay (ms)')
ylabel('Channels')
title('Planar Gradient Amplitude')
colorbar

%% Permutation test for signifiance with the method of maximum statistics
Nperm = 100;

% to account for the autocorrelation in both signals we implement the
% permutation with blocks of 10s
blocklen = 10 * Fs; % block length in samples
Nblock = ceil(Nt / blocklen);

% seperate the copula-normalised data into 10s blocks
bplnsum = cell(1,Nblock);
bspeech = cell(1,Nblock);
blen = zeros(Nblock,1);
for bi=1:Nblock
    idx = block_index(bi,blocklen,Nt);
    bplnsum{bi} = cplnsum(idx,:);
    bspeech{bi} = cspeech(idx);
    blen(bi) = length(idx);
end

% fix a delay to make this example quicker (would usually repeat the
% calcualtion over all delays for each permutation)
di = 5;
I = Iplnsum(:,di);
d = delay_samples(di);

Iperm = zeros(Nch,Nperm);
for pi=1:Nperm
    thsperm = randperm(Nblock);
    % apply delay/lag shift to shuffled blocks
    [dmeg, dspeech] = block_delay(bplnsum, bspeech(thsperm), d);
    for chi=1:Nch
        Iperm(chi, pi) = mi_gg(dmeg(:,chi), dspeech, true, true);
    end
end

% method of maximum statistics
% maximum values across permutations
Imax = max(Iperm,[],1);
thresh = prctile(Imax, 99);
Isig = I>thresh;

% if you have eeglab we can plot the topologies
if exist('topoplot')
    figure
    subplot(1,2,1)
    topoplot(I, chanlocs, 'maplimits', [0 max(I)]);
    title(sprintf('MI, lag = %d ms', delays(di)))
    colorbar
    
    subplot(1,2,2)
    topoplot(Isig, chanlocs, 'maplimits', [0 1]);
    title(sprintf('MI, p<0.01, lag = %d ms', delays(di)))
    colorbar
    
    colormap parula    
end

%% Multivariate GCMI: 2d vector

% We can calculate GCMI with multi-dimensional variables. Here we use the
% full 2d directional field vector (H and V)

% concatenate H and V along a new axis to create the bivariate response
pln2d = cat(2, reshape(plnH,[Nt 1 Nch]), reshape(plnV, [Nt 1 Nch]));

% Gaussian-copula normalisation, each variable, channel independently
cpln2d = copnorm(pln2d);
cspeech = copnorm(fltspc);


Ipln2d = zeros(Nch,Ndel);
for di=1:Ndel
    d = delay_samples(di);
    % resize speech according to the lag/delay considered
    dspeech = cspeech(1:(end-d));
    for chi=1:Nch
        % information calculation exacly as before except we pass in a 2d 
        % variable for MEG ( Nt x 2 )
        Ipln2d(chi,di) = mi_gg(cpln2d((1+d):end,:,chi), dspeech, true, true);
    end
end

figure
imagesc(delays,[],Ipln2d)
xlabel('MEG-Stimulus Delay (ms)')
ylabel('Channels')
title('Planar Gradient 2d Vector')
colorbar

%% Multivariate GCMI: direction of vector

% In the case of a multivariate directional vector we can normalise by 
% the amplitude to determin the MI available in the direction alone
% Here this is the direction of a magnetic fields, but the same approach
% can be applied to complex spectral data (taking real and imaginary part 
% as separate variables) to calculate the information conveyed in phase vs 
% the absolute value (power).

% concatenate H and V along a new axis to create the bivariate response
pln2d = cat(2, reshape(plnH,[Nt 1 Nch]), reshape(plnV, [Nt 1 Nch]));
% sum of square of compoments (vector amplitude)
plnsum = sqrt(plnH.^2 + plnV.^2);
% normalise away amplitude so points lie on unit circle (direction only)
plndir = bsxfun(@rdivide, pln2d, reshape(plnsum, [Nt 1 Nch]));

% Gaussian-copula normalisation, each variable, channel independently
cplndir = copnorm(plndir);
cspeech = copnorm(fltspc);


Iplndir = zeros(Nch,Ndel);
for di=1:Ndel
    d = delay_samples(di);
    % resize speech according to the lag/delay considered
    dspeech = cspeech(1:(end-d));
    for chi=1:Nch
        % information calculation exacly as before except we pass in a 2d 
        % variable for MEG ( Nt x 2 )
        Iplndir(chi,di) = mi_gg(cplndir((1+d):end,:,chi), dspeech, true, true);
    end
end

figure
imagesc(delays,[],Iplndir)
xlabel('MEG-Stimulus Delay (ms)')
ylabel('Channels')
title('Planar Gradient Direction')
colorbar
    