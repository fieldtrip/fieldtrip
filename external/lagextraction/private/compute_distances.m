function distances = compute_distances(points,mode)

% $Id: compute_distances.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

% disp(['    - Distance type = ',mode]);

npoints = size(points,1);

switch mode
    case 'eucli'
        % L_2 norm
        gram = points*points';
        norms = diag(gram);
        N = length(norms);
        distances = sqrt(repmat(norms,1,N) + repmat(norms',N,1) - 2*gram) ./ 2;
    case 'corr'
        % Correlation
        distances = sqrt((1 - corrcoef(points')) / 2);
    case 'corr2'
        % Squared Correlation
        distances = (1 - corrcoef(points')) / 2;
    case 'max'
        % distance of maximums
        maximums = max(points');
        distances = abs(repmat(maximums,npoints,1) - repmat(maximums',1,npoints));
    case 'inf'
        % norm inf
        distances = zeros(npoints);
        for i=1:N
            for j=(i+1):N
                d = norm(points(i,:)-points(j,:),'inf');
                distances(i,j) = size(points,1);
                distances(j,i) = distances(i,j);
            end
        end
    case 'histo'
        % L2 distance of histograms
        nb_bins = 20;
        minimum = min(points(:));
        maximum = max(points(:));
        bins = linspace(minimum,maximum,nb_bins);
        histograms = zeros(npoints,nb_bins);
        for i=1:npoints
            histograms(i,:) = hist(trial(i,:),bins);
        end
        gram = histograms*histograms';
        norms = diag(gram);
        distances = sqrt(repmat(norms,1,N) + repmat(norms',N,1) - 2*gram);
    case 'dgt'
        % L_2 norm of Discrete Gabor Transform
        % dgt_options.fmax = 50;
        dgt_options.fs = 256;
        dgt_options.doeven = 1;
        dgt_options.dt = 1;
        dgt_options.show = 0;
        dgt_options.tfr = 1;
        dgt_options.show = 0;
        [dcc,stfts] = perform_sgram(points',dgt_options);

        points = reshape(stfts,[],size(points,1))';
        gram = points*points';
        norms = abs(diag(gram));
        distances = sqrt(real(repmat(norms,1,npoints) + repmat(norms',npoints,1) - 2*gram));
    % case 'dtw'
    %     % dynamic time warpping
    %     distances = zeros(npoints);
    %     for k=1:nb_trials
    %         trial = points;
    %         [N,M] = size(trial);
    %         wrap_max = varargin{2};
    %         for i=1:N
    %             waitbar(k*i/(nb_trials*N),hw);
    %             for j=(i+1):N
    %                 d = dtwfast(trial(i,:),trial(j,:),wrap_max);
    %                 distances(i,j) = distances(i,j) + d;
    %                 distances(j,i) = distances(j,i) + d;
    %             end
    %         end
    %     end
    % case 'pdtw'
    %     % dynamic time warpping
    %     distances = zeros(npoints);
    %     for k=1:nb_trials
    %         trial = points;
    %         [N,M] = size(trial);
    %         wrap_max = varargin{2};
    %         for i=1:N
    %             waitbar(k*i/(nb_trials*N),hw);
    %             for j=(i+1):N
    %                 [dist,m,p,q,phi] = dtw(trial(i,:),trial(j,:),wrap_max);
    %                 d = sum(diag(m));
    %                 distances(i,j) = distances(i,j) + d;
    %                 distances(j,i) = distances(j,i) + d;
    %             end
    %         end
    %     end
    % case 'power'
    %     % spectral power or amplitude of wavelet coefficients
    %     N = npoints;
    %     distances = zeros(N);
    %     for i=1:N
    %         waitbar(i/N,hw);
    %         awts1 = cell(nb_trials);
    %         for k=1:nb_trials
    %             trial = points;
    %             awts1{k} = AWT(trial(i,:));
    %         end
    %         for j=(i+1):N
    %             awts2 = cell(nb_trials);
    %             for k=1:nb_trials
    %                 trial = points;
    %                 awts2{k} = AWT(trial(j,:));
    %             end
    %
    %             nb_scales = 15;
    %             nb_scales = 5;
    %             scales = (size(awts1{k},2)-15):size(awts1{k},2); % selecting a few low pass coefficients
    %             awts1{k} = awts1{k}(:,scales);
    %             awts2{k} = awts2{k}(:,scales);
    %             % ImageAWT(awts2{k},'Individual','jet');
    %
    %             d = 0;
    %             for k=1:nb_trials
    %                 d = d + norm(abs(awts1{k}).^2 - abs(awts2{k}).^2);
    %             end
    %             d = d / ( nb_trials * prod(size(awts1{k})) );
    %             distances(i,j) = d;
    %             distances(j,i) = d;
    %         end
    %     end
    %     for k=1:nb_trials
    %         trial = points;
    %         for i=1:N
    %             waitbar(k*i/(nb_trials*N),hw);
    %             for j=(i+1):N
    %                 awt1 = AWT(trial(i,:));
    %                 awt2 = AWT(trial(j,:));
    %                 d = norm(abs(awt1).^2-abs(awt2).^2);
    %                 distances(i,j) = distances(i,j) + d;
    %                 distances(j,i) = distances(j,i) + d;
    %             end
    %         end
    %     end
    % case 'phase'
    %     % phase locking using wavelet coeffecients
    %     N = npoints;
    %     distances = zeros(N);
    %     for i=1:N
    %         waitbar(i/N,hw);
    %         awts1 = cell(nb_trials);
    %         for k=1:nb_trials
    %             trial = points;
    %             awts1{k} = AWT(trial(i,:));
    %         end
    %         for j=(i+1):N
    %             awts2 = cell(nb_trials);
    %             for k=1:nb_trials
    %                 trial = points;
    %                 awts2{k} = AWT(trial(j,:));
    %             end
    %             d = 0;
    %             for k=1:nb_trials
    %                 tmp = awts1{k} ./ abs(awts1{k}) .* conj( awts2{k} ./ abs(awts2{k}) );
    %                 d = d + sum(tmp(:));
    %             end
    %             d = abs(d) / ( nb_trials * prod(size(awts1{k})) ); % d is the phase locking factor in [0,1]
    %             d = 1 - d; % set d = 0 when time series are similar
    %             distances(i,j) = d;
    %             distances(j,i) = d;
    %         end
    %     end
    otherwise
        error('Unknown distance type');
end

% close(hw)
