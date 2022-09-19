% This script provides a tutorial application of the Gaussian-Copula Mutual
% Information (GCMI) estimator with a 2 category discrete stimulus (face 
% vs noise images) event-related EEG data set.

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

% this will attempt to download the eeg data from the internet (~210MB)
% alternatively you can download the file eeg_face_vs_noise.mat manually
% and place it in the same directory as this script
fname = 'eeg_face_vs_noise.mat';
data_url = ['https://www.robince.net/data/gcmi/data/' fname];

% checks for data in current working directory
if ~exist(fname)
    disp(sprintf('Downloading %s ...', data_url))
    websave(fname,data_url);
    disp('Done.')
end

load(fname);
% csddat : [channels x time x trials]
% permute to trials first axis for looping
csddat = permute(csddat, [3 1 2]);
[Ntrl, Nch, Nt] = size(csddat);

% stim : stimulus class on each trial (0, 1)

%% Calculate GCMI at a specific sensor and time point

% this is the strongest effect
chi = 41;
ti = 168;

% for a t-test we contrast the data from the two classes
faceeeg = csddat(stim==0,chi,ti);
noiseeeg = csddat(stim==1,chi,ti);
[h,p,ci,stats] = ttest2(faceeeg, noiseeeg, 'VarType', 'unequal');
t = stats.tstat

% for MI we use the stimulus labels for each trial
% this works in the same way if there are more than 2 classes
% Reminder: stim must take values 0 or 1
I = gcmi_model_cd(csddat(:,chi,ti), stim, 2)


%% Calculate GCMI across all sensors and time point

% following the commonly used mass-univeriate approach we repeat 
% the calculation above for each sensor and time point
% we could call gcmi_cd as above inside the loop, but here we first
% normalise the data separately so we can reuse it for the permutation
% testing

% Gaussian-copula normalisation, works along first axis, applied to each 
% other dimension independently
ceeg = copnorm(csddat);
Ieeg = zeros(Nch,Nt);
for ti=1:Nt
    for chi=1:Nch
        % as the data has been copula-normalised we can use the 
        % Gaussian parametric estimator (whatever the distribution was
        % originally)
        Ieeg(chi,ti) = mi_model_gd(ceeg(:,chi,ti), stim, 2, true, true);
    end
end

figure
imagesc(time,[],Ieeg)
xlabel('Time')
ylabel('Channels')

%% Permutation test for signifiance with the method of maximum statistics
Nperm = 100;

% to reduce computation time, we consider a single time point here - the
% same method should normally be applied over all sensors and time points
ti = 168;
I = Ieeg(:,ti);
Iperm = zeros(Nch,Nperm);
for pi=1:Nperm
    idx = randperm(Ntrl);
    % randomly permute stimulus labels
    pstim = stim(idx);
    % repeat mass-univariate MI calculation
    for chi=1:Nch
        Iperm(chi,pi) = mi_model_gd(ceeg(:,chi,ti), pstim, 2, true, true);
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
    title(sprintf('MI, t = %d ms', time(ti)))
    colorbar
    
    subplot(1,2,2)
    topoplot(Isig, chanlocs, 'maplimits', [0 1])
    title(sprintf('MI, p<0.01, t = %d ms', time(ti)));
    colorbar
    
    colormap parula    
end



