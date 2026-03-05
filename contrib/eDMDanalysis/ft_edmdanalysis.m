function dataout = ft_edmdanalysis(cfg, datain)
%FT_EDMDANALYSIS performs eDMD decomposition on time series data
% 
% Copyright (c) 2026, David Chavez-Huerta
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any 
% later version.
%
%Use as:
%
%    [freq] = ft_edmdanalysis(cfg, raw)
%
%    OR (depending output options)
%
%    [raw] = ft_edmdanalysis(cfg, raw)
%
%The input data should be organised in a raw structure (like obtained from
%the FT_PREPROCESSING function). The output depends on the desired
%computation: either a reconstruction or a mode-based frequency analysis. 
%
%The configuration should contain:
%    cfg.output = different output options from the eDMD algorithm:
%                'raw': generates a reconstruction of the original data,
%                based on the mode decomposition 'seen through' the
%                dictionary and the first snapshot of the data.
%                'freq': generates a spectrum derived from Koopman 
%                eigenvalues
%                'binned_peak_freq': generates a freq structure with binned 
%                dominant frequencies. The bins can be customized and as 
%                default reflect brain rhythms. 
%    cfg        -  Configuration struct (see defaults below)
%   ---------------------------------------------------------------------
%    Optional inputs: 
%    cfg.cut        = double, rank truncation threshold control what 
%                     percentage of the total eigenvalue energy is preserved
%                     (default = 0.9999999999)
%    cfg.gamma      = double, RFF scaling parameter (default = 4)
%    cfg.D          = double, number of random Fourier features in the
%                     dictionary (default = 900)
%    cfg.nstacks    = double, Hankel stacking depth (default = 5)
%    cfg.MA         = double, moving average window (default = 0)
%    cfg.seed       = double, RNG seed for RFF reproducibility (default = 1)
%    cfg.freqEdges  = double vector, frequency bin edges. Default reflect 
%                     standard brain rhythms (default = [0 4 8 12 15 30 100])
%    cfg.verbose    = logical, print progress (default = true)
%    cfg.dictionary = char, {'rff','poly','hermite','identity'} dictionary 
%                     selection for the eDMD algorithm (default 'rff') 
%    cfg.poly_degree= double, used only if dictionary='poly' and sets the
%                     max polynomial degree. Increase carefully (default= 3)
%    cfg.hermite_degree = double, used only if dictionary='hermite' and sets
%                     the max degree. Increase carefully (default = 3)
%    cfg.foi          vector 1 x numfoi, frequencies of interest
%    cfg.smooth       0 or 1, optional smoothing in the frequency interpolation
%    cfg.normalize_recon = true/false (default = false) normalizes the
%                     reconstruction
%
%The configuration can also contain any parameters needed for the specific
%analysis:
%
%   cfg.channel
%   cfg.trials
%   cfg.latency 
%   etc...
%
%   ---------------------------------------------------------------------
%   ---------------------------------------------------------------------
%   Outputs:
%
%   The output structure depends on cfg.output:
%
%   ---------------------------------------------------------------------
%   cfg.output = 'freq'
%   ---------------------------------------------------------------------
%   Returns a FieldTrip-like frequency structure containing the interpolated
%   power spectrum derived from eDMD modes.
%
%   dataout.label        = {'edmd'} (virtual channel)
%   dataout.freq         = vector of frequencies (cfg.foi)
%   dataout.powspctrm    = [Ntrials x 1 x Nfreq] power spectrum
%   dataout.dimord       = 'rpt_chan_freq'
%   dataout.cfg          = configuration structure
%
%   Additional diagnostic outputs:
%   dataout.modefreqs    = cell array (1 x Ntrials) of mode frequencies
%   dataout.modepowers   = cell array (1 x Ntrials) of mode powers
%
%
%   ---------------------------------------------------------------------
%   cfg.output = 'binned_peak_freq'
%   ---------------------------------------------------------------------
%   Returns a FieldTrip-like frequency structure where power is aggregated
%   within predefined frequency bins (cfg.freqEdges). Peak frequencies per
%   bin are stored separately.
%
%   dataout.label        = {'edmd'} (virtual channel)
%   dataout.freq         = vector of bin centers
%   dataout.powspctrm    = [Ntrials x 1 x Nbins] summed power per bin
%   dataout.dimord       = 'rpt_chan_freq'
%   dataout.cfg          = configuration structure
%
%   Additional outputs:
%   dataout.peakfreq     = cell array (1 x Ntrials) of peak frequency per bin
%
%
%   ---------------------------------------------------------------------
%   cfg.output = 'raw'
%   ---------------------------------------------------------------------
%   Returns a FieldTrip raw structure containing the reconstructed signal
%   based on the Koopman mode decomposition.
%
%   dataout.label        = original channel labels
%   dataout.trial        = cell array (1 x Ntrials) with reconstructed data
%   dataout.time         = original time vectors
%   dataout.fsample      = sampling frequency
%   dataout.dimord       = 'rpt_chan_time'
%   dataout.cfg          = configuration structure
%
%
%   ---------------------------------------------------------------------
%   Notes:
%   ---------------------------------------------------------------------
%   - The 'freq' output corresponds to a spectrum obtained by
%     interpolating discrete Koopman mode frequencies.
%
%   - The 'binned_peak_freq' output provides a compact representation of
%     dominant oscillatory components within predefined frequency bands.
%
%   - The 'raw' output reconstructs the signal from the learned Koopman
%     model. Reconstruction is only performed when cfg.output = 'raw'.
%
%   - Mode-level information (modefreqs and modepowers) is always computed
%     internally but only returned for the 'freq' output.
% -----------------------------
% Preamble setup
% -----------------------------
ft_revision = '$Id$'; 
ft_nargin = nargin; 
ft_nargout = nargout; 
ft_preamble init 
ft_preamble debug 
ft_preamble loadvar datain 
ft_preamble provenance datain 
%ft_preamble_trackconfig          %removed from current FT version


% -----------------------------
% Validate input data
% -----------------------------
datain = ft_checkdata(datain, ...
  'datatype', {'raw','raw+comp'}, ...
  'hassampleinfo', 'yes');

cfg = ft_checkconfig(cfg, 'required', {});


% -----------------------------
% Set configuration defaults
% -----------------------------
cfg = ft_checkconfig(cfg, 'deprecated', {}); %Related to backwards compatibility
cfg = ft_checkconfig(cfg, 'renamed', {});    %Related to backwards compatibility
% Set defaults
cfg.cut            = ft_getopt(cfg, 'cut', 0.9999999999);
cfg.gamma          = ft_getopt(cfg, 'gamma', 4);
cfg.D              = ft_getopt(cfg, 'D', 900);
cfg.nstacks        = ft_getopt(cfg, 'nstacks', 5);
cfg.MA             = ft_getopt(cfg, 'MA', 0);
cfg.seed           = ft_getopt(cfg, 'seed', 1);
cfg.freqEdges      = ft_getopt(cfg, 'freqEdges', [0 4 8 12 15 30 100]);
cfg.verbose        = ft_getopt(cfg, 'verbose', true);
cfg.dictionary     = ft_getopt(cfg, 'dictionary', 'rff');   % 'rff' (default) 
cfg.poly_degree    = ft_getopt(cfg, 'poly_degree', 3);      % used only if dictionary='poly'
cfg.hermite_degree = ft_getopt(cfg, 'hermite_degree', 3);
cfg.foi            = ft_getopt(cfg, 'foi', [1:100]);
cfg.smooth         = ft_getopt(cfg, 'smooth', 0);
cfg.output         = ft_getopt(cfg, 'output', 'freq');
cfg.normalize_recon= ft_getopt(cfg, 'normalize_recon', false);

% Validate types of data and adequate intervals: 
cfg = ft_checkopt(cfg, 'dictionary', 'char', {'rff','poly','hermite','identity'});
cfg = ft_checkopt(cfg, 'foi', 'double');
cfg = ft_checkopt(cfg, 'output', 'char', {'freq','binned_peak_freq','raw'});
cfg = ft_checkopt(cfg, 'freqEdges', 'ascendingdoublevector');
if strcmp(cfg.dictionary, 'hermite') % validate hermite_degree only if hermite exists 
   cfg = ft_checkopt(cfg, 'hermite_degree', 'double');
end
if strcmp(cfg.dictionary, 'poly') % validate poly_degree only if poly exists 
   cfg = ft_checkopt(cfg, 'poly_degree', 'double');
end

if ~isnumeric(cfg.cut) || ~isscalar(cfg.cut) || cfg.cut<0 || cfg.cut>=1
   ft_error('cfg.cut must be numerical, scalar, inside the (0,1] interval');
end
if ~isnumeric(cfg.D) || ~isscalar(cfg.D) || cfg.D<=0
   ft_error('cfg.D must be scalar and positive');
end
if ~isnumeric(cfg.nstacks) ||   floor(cfg.nstacks)-cfg.nstacks ~= 0 || cfg.nstacks<1
   ft_error('cfg.nstacks must be numerical, no decimal part, equal or greater than 1');
end
if ~isnumeric(cfg.gamma) || ~isscalar(cfg.gamma) || cfg.gamma<1
   ft_error('cfg.gamma must be numerical, scalar, positive');
end

%Flag for reconstruction:
do_reconstruct = strcmp(cfg.output, 'raw');

% -----------------------------
% Channel / trial / latency selection
% -----------------------------
% Only use channel selection if field exists
if isfield(cfg,'channel')
  % Only keep channels that actually exist
   cfg.channel = intersect(cfg.channel, datain.label);
end

% we take only channel trials  and latency, so ft_selectdata doesn't flag
% the cfg.foi item:
cfgsel = [];   
cfgsel.channel = ft_getopt(cfg, 'channel', 'all');
cfgsel.trials  = ft_getopt(cfg, 'trials', 'all');
cfgsel.latency = ft_getopt(cfg, 'latency', 'all');

datain = ft_selectdata(cfgsel, datain);

% Validate again after selection
if isempty(datain.trial)
    ft_error('No trials left after selection. Check your cfg.channel, cfg.trials, and cfg.latency.');
end

% -----------------------------
% Define frequency grid
% -----------------------------

if isempty(cfg.foi)
    % default: 0 → Nyquist
    fmax = datain.fsample / 2;
    cfg.foi = linspace(0, fmax, 100);  % 100 points default
end

freq = cfg.foi(:)';   % ensure row freq vector
nFreq = length(freq); % freq cardinality


% -----------------------------
% Binning parameters
% -----------------------------

freqEdges = cfg.freqEdges;
numBins   = length(freqEdges)-1;

% -----------------------------
% Seed RNG for reproducibility
% -----------------------------

rng(cfg.seed, 'twister');

% -----------------------------
% Preallocate outputs
% -----------------------------

nTrials      = numel(datain.trial);
peakfeatures = cell(1, nTrials);
sumfeatures  = cell(1, nTrials);
modefreqs    = cell(1, nTrials);
modepowers   = cell(1, nTrials);
pow_all  = nan(nTrials, 1, nFreq);
xrecon = cell(1, nTrials);

% ----------------------
% Main eDMD loop: (process each trial independently)
% ----------------------
for tr = 1:nTrials
    X_zero = datain.trial{tr};   % channels x time (assumed, as in any FT function)
    if isempty(X_zero)
        if cfg.verbose, ft_warning('trial %d is empty -> skipping', tr); end
        continue
    end

    % X_zero dimensions
    [nChannels, nSamples] = size(X_zero);

    % optional moving average smoothing along time
    if cfg.MA > 0
        w = ones(1, cfg.MA)/cfg.MA;
        for ch = 1:nChannels
            X_zero(ch, :) = conv(X_zero(ch, :), w, 'same');
        end
    end

    % Build stacked data (Hankel-matrix)
    if cfg.nstacks > 1
        nst = cfg.nstacks;
        % create stacked matrix: (nChannels*nstacks) x newSamples
        newSamples = nSamples - nst + 1;
        Xaug = zeros(nChannels*nst, newSamples);
        for s = 1:nst
            Xaug((s-1)*nChannels + (1:nChannels), :) = X_zero(:, s:(s+newSamples-1));
        end
        X1 = Xaug(:, 1:end-1);
        X2 = Xaug(:, 2:end);
    else
        X1 = X_zero(:, 1:end-1);
        X2 = X_zero(:, 2:end);
    end

    % Update dims
    N = size(X1,1);   % dimension of stacked state
    M = size(X1,2);   % number of snapshots

    % ------------------------------------
    % Dictionary construction
    % ------------------------------------

    %Dictionary dimentional safeguards:

    %Error in case the user chooses an unreasonably large poly dictionary            
    if strcmp(cfg.dictionary,'poly')
       approx_terms = nchoosek(size(X1,1) + cfg.poly_degree, cfg.poly_degree);
    if approx_terms > 1e5
       ft_error('Polynomial dictionary is too large (>1e5 terms). Reduce degree or stacking.');
    end
    %Error in case the user chooses an unreasonably large rff dictionary
    end
    if strcmp(cfg.dictionary,'rff')
       feature_dim = 2 * cfg.D;
    if feature_dim > 5000
       ft_warning('Large RFF dictionary (%d features). Computation may be slow.', feature_dim);
    end
    end

    %Dictionary Switch: 
    switch lower(cfg.dictionary)

        case 'rff'
        D = cfg.D;
        W = randn(D, N) / cfg.gamma;
        b = 2*pi*rand(D,1);

        WX1 = W * X1 + b * ones(1, size(X1,2));
        WX2 = W * X2 + b * ones(1, size(X2,2));

        Psi_mX = (sqrt(1/D) * [cos(WX1); sin(WX1)])';
        Psi_mY = (sqrt(1/D) * [cos(WX2); sin(WX2)])';

        case 'poly'
        d = cfg.poly_degree;
        Psi_mX = buildPolynomialDictionary(X1, d);
        Psi_mY = buildPolynomialDictionary(X2, d);
        
        case 'hermite'
        d = cfg.hermite_degree;
        Psi_mX = buildHermiteDictionary(X1, d);
        Psi_mY = buildHermiteDictionary(X2, d);

        case 'identity' % equivalent to a DMD algorithm
        Psi_mX = X1';
        Psi_mY = X2';

        otherwise
        ft_error('Unknown dictionary type: %s', cfg.dictionary);
    end
    % ----------------------
    % Build compact K operator: Phi / B / U / S / V / Kt2 / 
    % ----------------------
      
    [U,S,V] = svd(Psi_mX / sqrt(M), 'econ');  
    singvals = diag(S).^2; %the sum of S square is the total energy
    r = find(cumsum(singvals) >= cfg.cut * sum(singvals), 1, 'first');
    rank_used(tr) = r;
    U = U(:,1:r);
    S = S(1:r,1:r);
    V = V(:,1:r);
    Kt2 = (S \ (U' * (Psi_mY/sqrt(M)) * V));    
    B = V*pinv(S)*U' * X1';   %using SVD of Psi_mX (faster matrix building)
    
    % ----------------------
    % Eigendecomposition of K: Wnn / MU / XI
    % ----------------------
    [XI,MU,W] = eig(Kt2);       %Generate eigenvalues and left-right eigenvectors
    inerp = sum(conj(W).*XI);                  %Calculate the w_k'xi_k products
    Wn = W./inerp;    %define the scaled w_n left eigenvectors in the u. circle
    scaledIP = sum(conj(Wn).*XI);              %Calculate the w_n'xi_k products
    ip_angle = angle(scaledIP);     %phase of w_n'xi_k products in the c. plane
    an_cor = exp(1i * ip_angle);                      %angle correction factor
    Wnn = Wn.*an_cor;                                 %normalized w eigevectors
    %ipfinal = sum(conj(Wnn).*XI); %this should be vector of ones if all went ok

    % ----------------------
    % Koompan modes, and state reconstruction: V_modes  / MU / XI
    % ----------------------

    %Calculate the Koopman modes v_i
    V_modes = (Wnn'*V'*B).';    %We must project l-eigenvectors on V' to recover time dimensions
    V_modes(isnan(V_modes)) = 0; %remove all undefined entries.
if do_reconstruct == 1   
    mu = repmat(diag(MU)', M, 1);
    row_indices = (1:M)';  % Create a column vector of row indices
    mu_p = mu .^ row_indices;  % Element-wise exponentiation
    Psi0 = Psi_mX(1,:);
    %State reconstruction: build the diagonal Phi matrix
    Phid = diag(Psi0*V*XI);  %We must project r-eigenvectors on V to recover time dimensions
    Xrr = real((mu_p*Phid*V_modes')'); %State recontruction: 
    orN = nChannels; %Original number of channels nChannels
    Xr = Xrr(1:orN,:); %Remove copies from the stacking 
    Xr_col = Xr(:);
if  cfg.normalize_recon
    Xr_col = normalize(Xr_col, "range");
end
    Xr_rec = reshape(Xr_col, orN, M);
end
    % ----------------------
    % Frequency and power calculation per mode
    % ----------------------
    % compute power per mode (norm squared):
    P_modes = vecnorm(V_modes).^2;
    % compute frequency per mode:
    %dt = 1; % user may change or pass in sampling freq if available (datain.fsample)
    dt = 1 / datain.fsample;
    % Frequencies derived from Koopman eigenvalues, they are not FFT-based
    % Absolute value used since spectrum is symmetric in ±f.
    f_modes = abs(angle(diag(MU)) / (2*pi*dt));
    f_modes(f_modes < 1e-6) = 0;     % threshold tiny freqs

    % Filter out zero-frequency modes 
    validIdx = f_modes > 0;
    fValid = f_modes(validIdx);
    PValid = P_modes(validIdx);

    % Sort & unique
    [fSorted, I] = sort(fValid);
    PSorted = PValid(I);
    [~, uniqueIdx] = unique(fSorted, 'stable');
    fUnique = fSorted(uniqueIdx);
    PUnique = PSorted(uniqueIdx);

    % Ensure column vectors
    fUnique = fUnique(:);
    PUnique = PUnique(:);

    % Remove NaNs or faulty entries
    valid = ~isnan(fUnique) & ~isnan(PUnique);
    fUnique = fUnique(valid);
    PUnique = PUnique(valid);
    
    % Interpolate frequencies into 'freq' array: 
if numel(fUnique) >= 2
    pow_interp = interp1(fUnique, PUnique, freq, 'linear', 0);
elseif numel(fUnique) == 1
    % single mode → assign nearest bin
    [~, idx] = min(abs(freq - fUnique));
    pow_interp = zeros(size(freq));
    pow_interp(idx) = PUnique;
else
    pow_interp = zeros(size(freq));
end

    % Optional smoothing 
if  isfield(cfg, 'smooth') && cfg.smooth > 0
    w = ones(1, cfg.smooth) / cfg.smooth;
    pow_interp = conv(pow_interp, w, 'same');
end

    pow_all(tr,1,:) = pow_interp;

    % Binning: compute peak frequency and sum power per bin
    peakvec = computeBinnedPeaks(fUnique, PUnique, freqEdges);
    sumvec  = computePowerSum(fUnique, PUnique, freqEdges);

    % Save outputs per trial
    peakfeatures{tr} = peakvec;
    sumfeatures{tr}  = sumvec;
    modefreqs{tr}    = fUnique;
    modepowers{tr}   = PUnique;
  if do_reconstruct == 1
    xrecon{tr}       = Xr_rec;
  end
   
    if cfg.verbose
    ft_info('trial %d processed ...', tr);
    end
end

%----------------------------------
% Wrap outputs into FT-style struct
%----------------------------------
nbins = length(freqEdges) - 1;
bin_centers = (freqEdges(1:end-1) + freqEdges(2:end)) / 2;

pow_bins = nan(nTrials, 1, nbins);

for tr = 1:nTrials
    if ~isempty(sumfeatures{tr})
        pow_bins(tr,1,:) = sumfeatures{tr};
    end
end

binnedout = [];

binnedout.label     = {'edmd'};
binnedout.freq      = bin_centers;
binnedout.powspctrm = pow_bins;
binnedout.dimord    = 'rpt_chan_freq';
binnedout.cfg = cfg;

% Keep peak frequencies as extra information
binnedout.peakfreq = peakfeatures;


dataout = [];
dataout.cfg = cfg;
dataout.peakfeatures = peakfeatures;
dataout.sumfeatures = sumfeatures;
dataout.modefreqs = modefreqs;
dataout.modepowers = modepowers;
dataout.rank = rank_used;

freqout = [];
freqout.label     = {'edmd'};
freqout.freq      = freq;
freqout.powspctrm = pow_all;
freqout.dimord    = 'rpt_chan_freq';
freqout.cfg       = cfg;

% Diagnostics (keep!)
freqout.modefreqs  = modefreqs;
freqout.modepowers = modepowers;
freqout.rank = rank_used;

rawout = [];
rawout.label   = datain.label;
rawout.trial   = xrecon;
% rawout.time needs a small correction since stacking consumes nstacks 
% time points
rawout.time = cell(1, nTrials);
for tr = 1:nTrials 
  if ~isempty(xrecon{tr})
    t0 = datain.time{tr}(1);
    rawout.time{tr} = t0 + (0:size(xrecon{tr},2)-1) / datain.fsample;
  else
    rawout.time{tr} = []; % or skip the trial earlier
  end
end
rawout.fsample = datain.fsample;
rawout.dimord  = 'rpt_chan_time';
rawout.cfg     = cfg;


%%%%%%%%%%%%%%%%%%%%%%%
% SWITCH OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.output

    case 'freq'
        
    freqout = ft_datatype_freq(freqout);
    dataout = freqout;

    case 'binned_peak_freq'

    binnedout = ft_datatype_freq(binnedout);
    dataout = binnedout;

    case 'raw'
        rawout = ft_datatype_raw(rawout);
        dataout = rawout;

    otherwise
        ft_error('Unsupported cfg.output: %s', cfg.output);
end
% ----------------------
% Postamble
% ----------------------
ft_postamble debug
ft_postamble previous datain
ft_postamble provenance dataout
ft_postamble history dataout
ft_postamble savevar dataout

end
% ----------------------
% Helper functions
% ----------------------
function innerbinpeaks = computeBinnedPeaks(f, P, freqEdges)
    nb = length(freqEdges)-1;
    innerbinpeaks = NaN(1, nb);
    for b = 1:nb
        idx = f >= freqEdges(b) & f < freqEdges(b+1);
        localpower = P(idx);
        localfreq  = f(idx);
        if isempty(localpower)
            innerbinpeaks(b) = NaN;
        else
            [~, Didx] = max(localpower);
            innerbinpeaks(b) = localfreq(Didx);
        end
    end
end

function powersum = computePowerSum(f, P, freqEdges)
    nb = length(freqEdges)-1;
    powersum = NaN(1, nb);
    for b = 1:nb
        idx = f >= freqEdges(b) & f < freqEdges(b+1);
        localpower = P(idx);
        if isempty(localpower)
            powersum(b) = NaN;
        else
            powersum(b) = sum(localpower);
        end
    end
end


function Psi = buildPolynomialDictionary(X, degree)
% Builds polynomial dictionary up to specified degree.

    [N, M] = size(X);

    % Start with constant term
    Psi = ones(1, M);

    % Add monomials degree 1 to "degree" value 
    for d = 1:degree
        combos = nchoosek(repmat(1:N,1,d), d);
        combos = unique(sort(combos,2),'rows');

        for k = 1:size(combos,1)
            term = ones(1,M);
            for j = 1:d
                term = term .* X(combos(k,j), :);
            end
            Psi = [Psi; term];
        end
    end

    Psi = Psi';  % matrix orientation (features x M)
end


function Psi = buildHermiteDictionary(X, degree)
% Hermite dictionary (probabilists' definition of Hermite polynomials)

    [N, M] = size(X);
    Xn = normalize(X, 2);  % normalize each channel across time for stability
    Psi = ones(1, M);  % H0
    for i = 1:N

        x = Xn(i,:);
        H_prev = ones(1,M);   % H0
        H_curr = x;           % H1
        Psi = [Psi; H_curr];
        for n = 1:(degree-1)
            H_next = x .* H_curr - n * H_prev;
            Psi = [Psi; H_next];
            H_prev = H_curr;
            H_curr = H_next;
        end
    end
    Psi = Psi';  % features x M
end