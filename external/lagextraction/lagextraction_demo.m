%% Lag extraction demo script
% this script runs an EEGLAB plugin used for lag extraction on single trial data.
%
% $Id: lagextraction_demo.m 4 2009-08-15 21:10:35Z gramfort $

%% Load data

USE_SYNTH = true;
USE_SYNTH = false;

clear options
close all

if USE_SYNTH

    randn('seed',0);

    clear EEG

    EEG.pnts = 500;
    EEG.srate = 1000; % 0.5 sec. recording
    EEG.trials = 100;
    EEG.nbchan = 1;
    EEG.comments = '';
    EEG.chanlocs(1).labels = 'Cz';
    EEG.chanlocs(1).X = 0;
    EEG.chanlocs(1).Y = 0;
    EEG.chanlocs(1).Z = 0;
    EEG.filename = 'synth';
    EEG.filepath = cd;

    EEG.times = linspace(0,EEG.pnts / EEG.srate,EEG.pnts);

    snr = 1.5;
    clear generate_options
    generate_options.snr = snr;

    if 1
        [EEG,ref,latencies] = generate_synth_eeg_data(EEG,generate_options);
    else
        clear generate_options

        generate_options.std_latency = 0.03;
        generate_options.mean_frequency = 3;
        generate_options.mean_latency = 0.3;

        % generate_options.mean_latency = 0.2;
        % generate_options.std_latency = 0.0;
        % generate_options.std_strech = 0.1;
        % generate_options.std_strech = 0.2;

        [EEG,ref,latencies] = generate_synth(EEG,generate_options);
    end
    EEG.times = linspace(0,1000 * EEG.pnts / EEG.srate,EEG.pnts);
    disp(['SNR : ',num2str(snr)]);

    channel = 1;
    time_win = [EEG.times(1) EEG.times(end)]; % (ms)
    bad_trials = []; % set bad trials
    use_ica = false;

else

    load('data/oddball3-num1-512Hz-chan10.set','-mat');

    channel = 1;
    time_win = [150 500]; % (ms)
    bad_trials = []; % set bad trials
    use_ica = false;

    latencies = eeg_getepochevent( EEG, {'200'}, [], 'latency');
    [tmp,motor_order] = sort(latencies,'descend');

    data = squeeze(EEG.data(1,:,:))';

    if  0
        % Hack to see result on motor response
        options.order = motor_order;
    end

end

%% Set parameters

options.distance = 'corr2';
% options.distance = 'corr'; % with corr

options.verbose = true; % show figures
options.disp_log = false; % print log messages

options.show_pca = true;
options.show_pca = false;

options.sigma = 0.01:0.01:0.2;

if 0
    % one alpha : no crossvalidation
    options.alpha = [0.1];
else
    % multiple alphas : crossvalidation (longer)
    options.alpha = [0.001,0.01,0.1];
end

%% Run extraction

[EEG, com, order, lags, event_type, E_lags] = pop_extractlag( EEG , use_ica, channel, time_win, options);

%% View results

if ~USE_SYNTH
    if options.verbose
        figure,
        imagesc(data(motor_order,:))
        title('Reordering with motor response')
        figure,
        imagesc(data(order,:))
        title('Reordering with Laplacian')
    end
end

if USE_SYNTH & options.verbose
    smart_figure('evoked_potentials');
    hold on
    plot(EEG.times,ref,'g','LineWidth',3);
    hold off
end

smart_figure('erp_data_ordered_lags');
if USE_SYNTH
    savefig(['lags_synth_alpha_',num2str(options.alpha)],22,{'pdf'})
else
    savefig(['lags_oddball_alpha_',num2str(options.alpha)],22,{'pdf'})
end
