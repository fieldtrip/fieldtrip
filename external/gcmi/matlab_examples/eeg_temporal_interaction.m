% This script provides a tutorial application of the Gaussian-Copula Mutual
% Information (GCMI) estimator to quantify temporal interactions in an
% event-related EEG response modulated by a continuous stimulus feature
% 
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

% this will attempt to download the eeg data from the internet (~4MB)
% alternatively you can download the file eeg_eye_visibility.mat manually
% and place it in the same directory as this script
fname = 'eeg_eye_visibility.mat';
data_url = ['https://www.robince.net/data/gcmi/data/' fname];

% checks for data in current working directory
if ~exist(fname)
    disp(sprintf('Downloading %s ...', data_url))
    websave(fname,data_url);
    disp('Done.')
end

load(fname);
% csddat : [trials x time]
[Ntrl, Nt] = size(csddat);
% eyestm : [trials x 2] (left, right eye visibility)

%% Calculate MI and correlation time course

% copula transform for GCMI
ceeg = copnorm(csddat);
ceyestim = copnorm(eyestim);

% repeat calculation at each time point
rho = zeros(1,Nt);
I = zeros(1,Nt);
for ti=1:Nt
    % correlation between left eye stimulus and EEG data at each time point
    rho(ti) = corr(csddat(:,ti), eyestim(:,1),'type', 'spearman');
    % GCMI between left eye stimulus and EEG data at each time point
    I(ti) = mi_gg(ceeg(:,ti), ceyestim(:,1), true, true);
end

figure
hold all
plot(time,rho)
plot(time,I)
legend('Spearmans rho', 'GCMI')

%% Conditional mutual information (CMI)

% note that there is a correlation between the left and right eye stimulus
% feature:
rho_lr = corr(eyestim(:,1),eyestim(:,2),'type','spearman')
% or with GCMI:
I_lr = mi_gg(ceyestim(:,1),ceyestim(:,2),true,true)

% we can address this with CMI which lets us measure the relationship
% between the left eye visibility and EEG signal, conditioning out
% variation of the right eye
CMI = zeros(1,Nt);
for ti=1:Nt
    % cmi function is like mi but takes three copula normalised inputs
    CMI(ti) = cmi_ggg(ceeg(:,ti), ceyestim(:,1), ceyestim(:,2), true, true);
end

figure
hold all
plot(time,I)
plot(time,CMI)
legend('MI: I(EEG; L)', 'CMI: I(EEG; L | R)')
Ifull = I;

%% Cross-temporal interaction information

% Computing interaction information between EEG(t1) and EEG(t2) requires
% calculating the MI in the joint response.
% All GCMI functions take multi-dimensional variables so we can calculate
% this just by calling MI with the two variables concatenated (samples
% first axis)
% for two time points:
t1 = 1;
t2 = 2;
Ijoint = gcmi_cc( [csddat(:,t1) csddat(:,t2)], eyestim(:,1) );
% or using pre-computed copula normalization as we will in a loop
Ijoint = mi_gg( [ceeg(:,t1) ceeg(:,t2)], ceyestim(:,1), true, true);

% To save time we only consider time points between 80 and 250 ms
inttidx = (time>80) & (time<250);
inttime = time(inttidx);
Nintt = length(inttime);

ceeg = copnorm(csddat(:,inttidx));
% we already calculated the MI at each time point:
I = Ifull(inttidx);

% To calculate the interaction matrix we loop over all pairs of points. The
% interaction information is symmetric so we only have to do one half
% triangle
Iint = zeros(Nintt,Nintt);
for t1=1:Nintt
    for t2=(t1+1):Nintt
        Ijoint = mi_gg( [ceeg(:,t1) ceeg(:,t2)], ceyestim(:,1), true, true );
        Iint(t1,t2) = Ijoint - I(t1) - I(t2);
    end
end
% fill in the symmetric other half
Iint = Iint + Iint';

figure
imagesc(inttime,inttime,Iint)
colormap parula
colorbar
axis square
xlabel('Time (ms)')
ylabel('Time (ms)')
title('Interaction information (bits)')

%% Calculating MI in the EEG gradient

% The synergy in the temporal interaction information reveals that at
% certain times the single trial temporal derivative / gradient of the EEG
% signal is modulated by the stimulus feature.

% first we calculate the temporal derivative of each trial
csdgrad = zeros(size(csddat));
for trli=1:Ntrl
    csdgrad(trli,:) = gradient(csddat(trli,:));
end

ceeg = copnorm(csddat);
ceeggrad = copnorm(csdgrad);
Igrad = zeros(1,Nt);
I2d = zeros(1,Nt);
for ti=1:Nt
    % we can calculate the MI in this signal as before
    % we use CMI throughout from now on to eliminate the effect of the
    % correlated right eye features
    Igrad(ti) = cmi_ggg(ceeggrad(:,ti), ceyestim(:,1), ceyestim(:,2), true, true);
    % we can also calculate the MI in the bivariate joint response of raw
    % signal + gradient... we just concatenate the matlab variables
    I2d(ti) = cmi_ggg( [ceeg(:,ti) ceeggrad(:,ti)], ceyestim(:,1), ceyestim(:,2), true, true);
end

figure
hold all
plot(time,CMI)
plot(time,Igrad)
plot(time,I2d)
legend('Raw voltage', 'Gradient', '2D raw and gradient')
xlabel('Time (ms)')
ylabel('CMI (bits)')

%% Calculating interaction information in the 2d gradient response

% We can calculate cross-temporal interaction information exactly as before
% but using our 2d voltage + gradient response, and also conditioning out
% the confounding right eye stimulus

ceeg = copnorm(csddat(:,inttidx));
ceeggrad = copnorm(csdgrad(:,inttidx));

% we already calculated the 2D CMI at each time point:
I = I2d(inttidx);

% Interaction matrix
Iint2d = zeros(Nintt,Nintt);
for t1=1:Nintt
    for t2=(t1+1):Nintt
        % now we concatenate the 4 variables: voltage + gradient at each
        % time point
        Ijoint = cmi_ggg( [ceeg(:,t1) ceeggrad(:,t1) ceeg(:,t2) ceeggrad(:,t2)],...
                           ceyestim(:,1), ceyestim(:,2), true, true );
        Iint2d(t1,t2) = Ijoint - I(t1) - I(t2);
    end
end
% fill in the symmetric other half
Iint2d = Iint2d + Iint2d';

figure
% plot interaction information
subplot(1,2,1)
imagesc(inttime,inttime,Iint2d)
colormap parula
colorbar
cl = max(abs(caxis));
caxis([-cl cl])
axis square
xlabel('Time (ms)')
ylabel('Time (ms)')
title('Interaction information (bits)')
% plot redundancy only for more detailed view
subplot(1,2,2)
red2d = -Iint2d;
imagesc(inttime,inttime,red2d);
colormap parula
colorbar
cl = max(abs(caxis));
caxis([0 cl])
axis square
xlabel('Time (ms)')
ylabel('Time (ms)')
title('Redundancy (bits)')

%% Calculating the emergence of novel MI over time

% We downsample the data for this analysis:
% downsample the raw signal
dscsd = resample(csddat',1,4)';
% downsample the gradient (calculated with full temporal resolution)
dscsdgrad = resample(csdgrad',1,4)';
dstime = resample(time,1,4);
dsNt = length(dstime);

ceeg = copnorm(dscsd);
ceeggrad = copnorm(dscsdgrad);

% calculate CMI as before for the downsampled data
I2d = zeros(1,dsNt);
for ti=1:dsNt
    I2d(ti) = cmi_ggg( [ceeg(:,ti) ceeggrad(:,ti)], ceyestim(:,1), ...
                        ceyestim(:,2), true, true);
end

% calculate new MI at each time point
I2dnew = zeros(1,dsNt);
for ti=2:dsNt
    % add extra variables to condition out with Matlab concatenation
    % as well as the right eye we also condition out the 2d response at 
    % the previous time point
    I2dnew(ti) = cmi_ggg( [ceeg(:,ti) ceeggrad(:,ti)], ceyestim(:,1), ...
                          [ceeg(:,ti-1) ceeggrad(:,ti-1) ceyestim(:,2)], true, true);
end

% manual inspection of the I2dnew time course reveals two clear early peaks
% to ensure later novel information is genuinely new we include these peaks
% into the conditional
peaks = [36 42]; % manually identified from I2dnew
I2dnewpeak = zeros(1,dsNt);
for ti=2:dsNt
    % construct multi-dimensional variable to condition out
    condvar = [ceeg(:,ti-1) ceeggrad(:,ti-1) ceyestim(:,2)];
    for pi=1:length(peaks)
        if ti>(peaks(pi)+1)
            % if this peak is earlier than the current time point
            % add it to the conditioning
            condvar = [condvar ceeg(:,peaks(pi)) ceeggrad(:,peaks(pi))];
        end
    end
    I2dnewpeak(ti) = cmi_ggg( [ceeg(:,ti) ceeggrad(:,ti)], ceyestim(:,1), ...
                               condvar, true, true);
end

figure
hold all
plot(dstime,I2d)
plot(dstime,I2dnew)
plot(dstime,I2dnewpeak)
xlabel('Time (ms)')
ylabel('bits')
legend('MI','New MI: previous timepoint', 'New MI: previous timepoint + peaks')
