function test_ft_edmdanalysis_md
% Ensure FieldTrip defaults are here
ft_defaults

fprintf('Running test_ft_edmdanalysis...\n')

%% ------------------------------------------------------------------------
% Create synthetic MULTI-CHANNEL oscillatory dataset
% -------------------------------------------------------------------------

fs = 500;                 
T  = 2;                   
t  = 0:1/fs:T-1/fs;       
f0 = 10;                  

rng(1)  % deterministic noise

signal = sin(2*pi*f0*t);

nchan = 4;
X = zeros(nchan, length(t));

X(1,:) = signal;                 % oscillation in channel 1
X(2,:) = 0.05 * randn(size(t));  % small noise
X(3,:) = 0.05 * randn(size(t));
X(4,:) = 0.05 * randn(size(t));

data = [];
data.label      = {'chan1','chan2','chan3','chan4'};
data.fsample    = fs;
data.trial      = {X};
data.time       = {t};
data.sampleinfo = [1 numel(t)];


% ------------------------------------------------------------------------
%1) TEST FREQ OUTPUT
%-------------------------------------------------------------------------

cfg = [];
cfg.output     = 'freq';
cfg.dictionary = 'hermite';   % equivalent to DMD
cfg.nstacks    = 4;            % no Hankel stacking
cfg.foi        = 0:0.5:40;     % frequency grid
cfg.verbose    = false;

freq = ft_edmdanalysis(cfg, data);

% ---- Basic structural checks
assert(isfield(freq, 'powspctrm'), 'Missing powspctrm field')
assert(strcmp(freq.dimord, 'rpt_chan_freq'), 'Wrong dimord')
assert(size(freq.powspctrm,1) == 1, 'Wrong number of trials')
assert(size(freq.powspctrm,2) == 1, 'Wrong number of channels')

% ---- Check frequency detection
[~, idx] = max(squeeze(freq.powspctrm));
f_detected = freq.freq(idx);

assert(abs(f_detected - f0) < 1.0, ...
    'Detected frequency deviates too much from ground truth.')

fprintf('✓ freq output test passed (detected %.2f Hz)\n', f_detected)


% ------------------------------------------------------------------------
% 2) TEST BINNED OUTPUT
% -------------------------------------------------------------------------

cfg = [];
cfg.output     = 'binned_peak_freq';
cfg.dictionary = 'hermite';
cfg.nstacks    = 4;
cfg.freqEdges  = [0 4 8 12 20 40];
cfg.verbose    = false;

binned = ft_edmdanalysis(cfg, data);

assert(isfield(binned, 'powspctrm'), 'Missing binned powspctrm')
assert(strcmp(binned.dimord, 'rpt_chan_freq'), 'Wrong dimord (binned)')
assert(length(binned.freq) == length(cfg.freqEdges)-1, ...
    'Wrong number of bins')

% 10 Hz should fall in 8–12 Hz bin
[~, maxbin] = max(squeeze(binned.powspctrm));
bin_center = binned.freq(maxbin);

assert(bin_center >= 8 && bin_center <= 12, ...
    'Peak not detected in correct frequency bin.')

fprintf('✓ binned_peak_freq test passed (bin center %.2f Hz)\n', bin_center)


% ------------------------------------------------------------------------
% 3) TEST RAW RECONSTRUCTION OUTPUT
% -------------------------------------------------------------------------

cfg = [];
cfg.output     = 'raw';
cfg.dictionary = 'hermite';
cfg.nstacks    = 4;
cfg.verbose    = false;

rawrec = ft_edmdanalysis(cfg, data);

assert(isfield(rawrec, 'trial'), 'Missing trial field in raw output')
assert(strcmp(rawrec.dimord, 'rpt_chan_time'), 'Wrong raw dimord')

% reconstruction length should be original-1 due to X1/X2 split
expected_length = length(t) - cfg.nstacks;
rec_length = size(rawrec.trial{1},2);

assert(rec_length == expected_length, ...
    'Unexpected reconstruction length.')

fprintf('✓ raw reconstruction test passed (length %d)\n', rec_length)


%% ------------------------------------------------------------------------
fprintf('All ft_edmdanalysis tests PASSED.\n')

end