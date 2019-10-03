% pop_extractlag() - If less than
%           4 arguments are given, a window pops up
%           to ask for the value of the additional
%           parameters.
%
% Usage:
%    >>    [OUTEEG, COM, ORDER, LAGS_MS, EVENT_TYPE] = pop_extractlag( INEEG, USE_ICA, CHANNEL_ID, TIME_WIN, OPTIONS );
%
% Inputs:
%    INEEG      - input EEG dataset
%    use_ica    - type of processing. 1 process the raw data and 0 the ICA components.
%    channel_id - channel or component id to consider.
%    time_win   - time window in ms on which lag extraction should be performed
%    options    : structure containing options with default values
%        Some important options :
%           options.use_maximum = true (resp. false) : lag pass through a maximum (resp. lag pass through a minimum)
%           options.sigma = 0.05 : controls the similarity measures between time series
%           options.alpha = 0.01 (example of value) : controls regularity of lag function
%
% Outputs:
%    OUTEEG    - output dataset
%
% See also:
%    EEGLAB

% (c) Copyright 2008 Alexandre Gramfort INRIA. All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA    02111-1307    USA

% $Id: pop_extractlag.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

function [EEG, com, order, lags_ms, event_type, E_lags] = pop_extractlag( EEG, use_ica, channel_id, time_win, options )

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
    % if the user press the cancel button

% display help if not enough arguments
% ------------------------------------
if nargin < 1
    help pop_extractlag;
    return;
end;

% Set default options

if nargin<5
    options.null = 0;
end

if ~isfield(options, 'NORMALIZE_DATA')
    options.NORMALIZE_DATA = 1;
end
NORMALIZE_DATA = options.NORMALIZE_DATA;

if ~isfield(options, 'CENTER_DATA')
    options.CENTER_DATA = 0;
end
CENTER_DATA = options.CENTER_DATA;

if ~isfield(options, 'USE_ADAPTIVE_SIGMA')
    options.USE_ADAPTIVE_SIGMA = 0;
end
USE_ADAPTIVE_SIGMA = options.USE_ADAPTIVE_SIGMA;

if ~isfield(options, 'DISPLAY_SIGNAL')
    options.DISPLAY_SIGNAL = 0;
end
DISPLAY_SIGNAL = options.DISPLAY_SIGNAL;

if ~isfield(options, 'sigma')
    sigma_txt = '0.01:0.01:0.2';
    options.sigma = eval([ '[' sigma_txt ']' ]);
end
sigma = options.sigma;

if ~isfield(options, 'distance')
    options.distance = 'corr2';
end
distance = options.distance;

if ~isfield(options, 'CONNECT_EMBEDDED_POINTS')
    options.CONNECT_EMBEDDED_POINTS = 0;
end
CONNECT_EMBEDDED_POINTS = options.CONNECT_EMBEDDED_POINTS;

if ~isfield(options, 'alpha')
    alpha_txt = '0.001,0.01,0.1,1';
    options.alpha = eval([ '[' alpha_txt ']' ]);
end
alpha = options.alpha;

if ~isfield(options, 'exponent')
    options.exponent = 1;
end
exponent = options.exponent;

if ~isfield(options, 'use_maximum')
    options.use_maximum = true;
end
use_maximum = options.use_maximum;

if ~isfield(options, 'show_pca')
    options.show_pca = false;
end
show_pca = options.show_pca;

if ~isfield(options, 'show_trial_number')
    options.show_trial_number = false;
end
show_trial_number = options.show_trial_number;

if ~isfield(options, 'verbose')
    options.verbose = true;
end
verbose = options.verbose;

if ~isfield(options, 'disp_log')
    options.disp_log = true;
end
disp_log = options.disp_log;

% pop up window
% -------------

if nargin < 4

    cb_plotcomps = [ 'pop_topoplot(EEG);' ];
    cb_erpimage = [ 'pop_erpimage(EEG,1);' ];

    listui = { { 'style' 'text'       'string' 'Use ICA components' } ...
               { 'style' 'checkbox'   'string' '' } ...
               { 'style' 'text'       'string' '(check=yes)' } ...
               { 'style' 'text'       'string' 'Channel or Component' } ...
               { 'style' 'edit'       'string' [ num2str(1) ] } ...
               { 'style' 'pushbutton' 'string' 'Plot topo.' 'callback' cb_plotcomps }  ...
               { 'style' 'text'       'string' 'Time window (ms)' } ...
               { 'style' 'edit'       'string' [ num2str(EEG.times(1)) ':' num2str(EEG.times(end)) ] } ...
               { 'style' 'pushbutton' 'string' 'Erpimage' 'callback' cb_erpimage } ...
               { 'style' 'text'       'string' 'Use minimum' } ...
               { 'style' 'checkbox'   'string' '' } ...
               { 'style' 'text'       'string' '(check=yes)' } ...
               { 'style' 'text'       'string' 'Sigma' } ...
               { 'style' 'edit'       'string' sigma_txt } ...
               { 'style' 'text'       'string' '' } ...
               { 'style' 'text'       'string' 'Alpha' } ...
               { 'style' 'edit'       'string' alpha_txt } ...
               { 'style' 'text'       'string' '' } ...
               { 'style' 'text'       'string' 'Verbose' } ...
               { 'style' 'checkbox'   'string' '' } ...
               { 'style' 'text'       'string' '(check=yes)' } ...
             };
    geom = { [1 1 0.5] [0.8 0.2 1] [1 1 1] [1 1 1] [1 1 1] [1 1 1] [1 1 1] };

    results = inputgui('geometry', geom, 'uilist', listui, 'helpcom', 'pophelp(''pop_extractlag'')', ...
        'title', 'Extract lags');

    if isempty(results), return; end;

    use_ica             = eval( [ '[' num2str(results{1}) ']' ] );
    channel_id          = eval( [ '[' results{2} ']' ] );
    time_win            = eval( [ '[' results{3} ']' ] );
    time_win = [time_win(1),time_win(end)];

    options.use_maximum = not(results{4});
    options.sigma       = eval( [ '[' results{5} ']' ] );
    options.alpha      = eval( [ '[' results{6} ']' ] );
    options.verbose     = results{7};

end;

tmin = find(EEG.times > time_win(1),1,'first');
tmax = find(EEG.times < time_win(2),1,'last');

% call function
% ---------------------------------------------------
if size(EEG.data,3) == 1
    error('EEG dataset is not epoched (EEG.data should of dimensions : nb channels x nb times x nb epochs)')
end

if use_ica == 0
    data = squeeze(EEG.data(channel_id,:,:))';
    points = squeeze(EEG.data(channel_id,tmin:tmax,:))';
else
    if ~isempty( EEG.icadata )
        data = squeeze(EEG.icadata(channel_id,:,:))';
        points = squeeze(EEG.icadata(channel_id,tmin:tmax,:))';
    else
        error('You must run ICA first');
    end;
end;
npoints = size(data,1);

%% Cross Validation
alphas = options.alpha;

if length(alphas) > 1 % Use Cross validation error if multiple alphas
    best_CVerr = -Inf;

    K = 10;
    disp(['--- Running K Cross Validation (K = ',num2str(K),')']);

    block_idx = fix(linspace(1,npoints,K+1)); % K cross validation
    for jj=1:length(alphas)
        options.alpha = alphas(jj);

        CVerr = 0;
        for kk=1:K
            bidx = block_idx(jj):block_idx(jj+1);
            idx = 1:npoints;
            idx(bidx) = [];

            data_k = data(idx,:);
            points_k = points(idx,:);
            [order,lags,coords,E_lags] = extractlag(points_k,options);

            data_reordered = data_k(order,:);
            lags = lags+tmin;
            [data_aligned,data_aligned_times] = perform_realign( data_reordered, EEG.times, lags );
            ep_evoked = mean(data_aligned);
            ep_evoked(isnan(ep_evoked)) = 0;
            ep_evoked(ep_evoked==Inf) = 0;

            ep_evoked = ep_evoked ./ norm(ep_evoked);

            data_k = data(bidx,:);
            data_k = normalize_rows(data_k);

            ep_raw = mean(data);
            ep_raw = ep_raw ./ norm(ep_raw);
            for pp=1:length(bidx)
                c = xcorr(ep_raw,data_k(pp,:));
                % max(c(:))
                c = xcorr(ep_evoked,data_k(pp,:));
                % max(c(:))
                [mc,mci] = max(c(:));
                % keyboard

                CVerr = CVerr + max(c(:));
            end
        end

        CVerr = CVerr/npoints;

        if CVerr > best_CVerr
            best_CVerr = CVerr;
            best_alpha = alphas(jj);
        end
    end

    options.alpha = best_alpha;
end

if use_maximum
    [order,lags,coords,E_lags] = extractlag( points, options );
else
    [order,lags,coords,E_lags] = extractlag( -points, options );
end

disp(['---------- Using alpha = ',num2str(options.alpha)]);

data_reordered = data(order,:);
lags = lags+tmin;
[data_aligned,data_aligned_times] = perform_realign( data_reordered, EEG.times, lags );

order_inv = perminv(order);
lags_no_order = lags(order_inv);
lags_ms = EEG.times(lags_no_order);

event_type = [num2str(time_win(1)),' - ',num2str(time_win(2)),' ms'];

if isfield(EEG,'chanlocs') % Real dataset (not synthetic)
    event_type = [EEG.chanlocs(channel_id).labels,' (',num2str(channel_id),') : ',event_type];
end

if use_maximum
    event_type = [event_type ' (max)'];
else
    event_type = [event_type ' (min)'];
end

% try

    clear events
    events = zeros(length(lags),3);

    for ii=1:length(lags)

        clear event
        event.type = event_type;
        event.latency = (ii-1)*length(EEG.times) + lags_no_order(ii);

        if ~isfield(EEG,'urevent')
            nurevent = 0;
            EEG.urevent(1) = event;
        else
            nurevent = length(EEG.urevent);
            EEG.urevent(end+1) = event;
        end

        event.code = 'Realignement';
        event.duration = 1;
        % event.channel = channel_id; % or 0 for all channels
        event.channel = 0; % or 0 for all channels
        event.time = [];

        % event.urevent = nurevent+1;
        event.epoch = ii;

        if isfield(EEG,'event')
            evidx = length(EEG.event)+1;
        else
            evidx = 1;
        end
        EEG.event(evidx).type = event.type;
        EEG.event(evidx).latency = event.latency;
        EEG.event(evidx).epoch = event.epoch;

        if ~isfield(EEG,'event')
            EEG.event(1) = event;
        else
            if isfield(EEG.event,'code'), EEG.event(evidx).code = event.code; end
            if isfield(EEG.event,'duration'), EEG.event(evidx).duration = event.duration; end
            if isfield(EEG.event,'channel'), EEG.event(evidx).channel = event.channel; end
            if isfield(EEG.event,'time'), EEG.event(evidx).time = event.time; end
            if isfield(EEG.event,'epoch'), EEG.event(evidx).epoch = event.epoch; end
        end
    end

    EEG = eeg_checkset(EEG);

% catch
%     warning('Problem with EEG set structure : Markers could not be correctly added to the dataset !')
% end

% return the string command
% -------------------------
com = sprintf('pop_extractlag( %s, %d, [%s], [%s] );', inputname(1), use_ica, int2str(channel_id), num2str(time_win));

%% view embedding results

t0 = find(EEG.times > 0,1,'first');

if options.verbose
    smart_figure('data');
    % plot(EEG.times,data(fix(linspace(1,size(data,1),4)),:));
    plot(EEG.times,data(fix(linspace(1,size(data,1),2)),:));
    axis tight
    xlabel('Time (ms)')
    line([time_win(1) time_win(1)],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([time_win(2) time_win(2)],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([0 0],get(gca,'ylim'),'Color','k','LineWidth',2,'LineStyle','-');

    smart_figure('erp_data');
    imagesc(data,[min(data(:)),max(data(:))]);
    % colorbar
    xlabel('Time (ms)')
    ylabel('Trial')
    title('Data')
    set(gca,'XTickLabel',fix(EEG.times(get(gca,'XTick'))))
    line([tmin tmin],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([tmax tmax],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([t0 t0],get(gca,'ylim'),'Color','k','LineWidth',2,'LineStyle','-');

    smart_figure('erp_data_ordered');
    imagesc(data_reordered,[min(data(:)),max(data(:))]);
    % colorbar
    xlabel('Time (ms)')
    ylabel('Trial')
    title('Data Ordered')
    set(gca,'XTickLabel',fix(EEG.times(get(gca,'XTick'))))
    % set(gcf,'WindowButtonDownFcn','callbackRaster(EEG)'); % set the callback
    line([tmin tmin],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([tmax tmax],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([t0 t0],get(gca,'ylim'),'Color','k','LineWidth',2,'LineStyle','-');

    smart_figure('2D_embedding');
    plot(coords(:,1),coords(:,2),'+');
    xlabel('y_1')
    ylabel('y_2')
    if show_trial_number
        dxtext = std(coords(:,1))*0.1;
        dytext = std(coords(:,2))*0.1;
        for i=1:size(coords,1)
            text(coords(i,1)+dxtext,coords(i,2)+dytext,num2str(i));
        end
    end
    set(get(gca,'Title'),'String','2D Embedding');
end

if options.show_pca
    % [coeff, SCORE, LATENT] = princomp(points(order,:));
    [coeff, SCORE, LATENT] = princomp(zscore(points(order,:)));

    smart_figure('PCA embedding'); clf
    scatter(SCORE(:,1),SCORE(:,2),25,jet(npoints),'filled');
    % scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),25,jet(npoints),'filled'); cameramenu
    axis square
end

%% view lag extraction results

if 1 & options.verbose
    smart_figure('erp_data_ordered_lags');
    imagesc(data_reordered,[min(data(:)),max(data(:))]);
    % colorbar
    xlabel('Time (ms)')
    ylabel('Trial')
    title('Data Ordered')
    set(gca,'XTickLabel',fix(EEG.times(get(gca,'XTick'))))
    % set(gcf,'WindowButtonDownFcn','callbackRaster(EEG)'); % set the callback
    line([tmin tmin],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([tmax tmax],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
    line([t0 t0],get(gca,'ylim'),'Color','k','LineWidth',2,'LineStyle','-');
    hold on % super impose the lags and the smoothed maximums
    plot(lags,[1:npoints],'k+')
    plot(lags,[1:npoints],'k')
    line([mean(lags) mean(lags)],get(gca,'ylim'),'Color','b','LineWidth',1,'LineStyle','--')
    hold off
end

if 1 & options.verbose
    smart_figure('erp_data_aligned');
    imagesc(data_aligned,[min(data_aligned(:)),max(data_aligned(:))])
    xlabel('Time (ms)')
    ylabel('Trial')
    title('Data Ordered and aligned')
    set(gca,'XTickLabel',fix(EEG.times(get(gca,'XTick'))))
    line([mean(lags) mean(lags)],get(gca,'ylim'),'Color','b','LineWidth',1,'LineStyle','--')
    line([t0 t0],get(gca,'ylim'),'Color','k','LineWidth',2,'LineStyle','-');
end

%% view evoked potential after realigning

ep_orig = mean(data_reordered); % evoked potential init
ep_aligned = mean(data_aligned); % evoked potential after realignment

smart_figure('evoked_potentials');
plot(EEG.times,ep_orig,'r');
hold on
plot(EEG.times,ep_aligned,'b');
hold off
xlabel('Time (ms)')
ylim = get(gca,'ylim');
ylim(1) = ylim(1) - 0.1 * abs(ylim(2)-ylim(1));
ylim(2) = ylim(2) + 0.1 * abs(ylim(2)-ylim(1));
set(gca,'ylim',ylim);
title('Naive averaging vs realigned averaging')
line([time_win(1) time_win(1)],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
line([time_win(2) time_win(2)],get(gca,'ylim'),'Color','k','LineWidth',1,'LineStyle','--');
line([0 0],get(gca,'ylim'),'Color','k','LineWidth',2,'LineStyle','-');
return;
