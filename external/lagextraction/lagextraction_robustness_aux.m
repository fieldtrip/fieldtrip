pref = 'gc_';

for kk=1:length(SNRs)
    randn('seed',rseed);
    rand('seed',rseed);
    % randn('seed',0);
    % randn('seed',sum(100*clock));

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
    EEG.xmax = 0.5;
    EEG.xmin = 0;
    EEG.times = linspace(0,EEG.pnts / EEG.srate,EEG.pnts);

    clear generate_options
    if USE_AR
        generate_options.use_ar = true;
    end

    switch template_id
    case 1
        generate_options.snr = SNRs(kk);
        EEG.times = linspace(0,EEG.pnts / EEG.srate,EEG.pnts);
        [EEG,ref,latencies] = generate_synth_eeg_data(EEG,generate_options);
        EEG.times = linspace(0,1000 * EEG.pnts / EEG.srate,EEG.pnts);
    case 2
        clear generate_options
        % generate_options.std_latency = 0.03;
        generate_options.std_latency = 0.05;
        generate_options.mean_frequency = 3;
        generate_options.mean_latency = 0.3;
        generate_options.snr = SNRs(kk);
        [EEG,ref,latencies] = generate_synth(EEG,generate_options);
    case 3
        clear generate_options
        % generate_options.std_latency = 0.03;
        generate_options.std_latency = 0.05;
        generate_options.mean_frequency = 15;
        generate_options.mean_latency = 0.3;
        generate_options.snr = SNRs(kk);
        [EEG,ref,latencies] = generate_synth(EEG,generate_options);
    end

    disp(['SNR : ',num2str(SNRs(kk))]);

    % lagextraction parameters
    channel = 1;
    time_win = [EEG.times(1) EEG.times(end)]; % (ms)
    bad_trials = []; % set bad trials
    use_ica = false;

    clear options
    options.distance = 'corr2';
    options.sigma = linspace(0.01,0.2,15);
    options.alpha = [0,0.001,0.01,0.1];
    options.verbose = false;
    [EEG, com, order, lags] = pop_extractlag( EEG , use_ica, channel, time_win, options);

    data_orig = squeeze(EEG.data(channel,:,:))';
    EEG = pop_epoch( EEG, {EEG.event(end).type}, [-0.15 0.15]);
    ep_aligned = mean(squeeze(EEG.data(channel,:,:)),2);
    times_aligned = EEG.times;

    ep_aligned_norm = ep_aligned ./ norm(ep_aligned);
    ref_norm = ref ./ norm(ref);
    [tmp,idx]  = max(xcorr(ep_aligned,ref));
    ref_norm = circshift(ref_norm(:),idx);
    ref_norm = ref_norm(1:length(ep_aligned_norm));
    cref(kk) = corr(ref_norm./norm(ref_norm),ep_aligned_norm);

    latencies = latencies - mean(latencies);
    data_gt = [];
    for pp=1:EEG.trials
        data_gt(pp,:) = circshift(data_orig(pp,:)',-fix(EEG.srate*latencies(pp)));
    end
    ep_gt = mean(data_gt);
    ep_gt = ep_gt ./ norm(ep_gt);

    [tmp,idx]  = max(xcorr(ep_aligned,ep_gt));
    ep_gt = circshift(ep_gt(:),idx);
    ep_gt = ep_gt(1:length(ep_aligned_norm));

    cref_gt(kk) = corr(ref_norm./norm(ref_norm),ep_gt./norm(ep_gt));

    % close all
    % plot(circshift(ref_norm(:),idx))
    % hold on
    % plot(ep_aligned_norm)
    % hold off

    % plot([ref_norm(:),ep_aligned_norm(:)])
    % return

    if options.verbose
        smart_figure('evoked_potentials');
        hold on
        plot(times_orig,ref,'g','LineWidth',3);
        hold off
        pause
    end
end
