function [fs, fg, opt] = ft_specest_fitparametric(powspctrm, freqs, opt, method)

if ~isequal(method, 'fooof')
  ft_error('methods other than fooof are not supported');
end

% check for brainstorm functions on the path, and add if needed
ft_hastoolbox('brainstorm', 1);

% This grabs the defaults from the brainstorm code
opts_bst  = getfield(process_fooof('GetDescription'), 'options');

% Added experimental functionality to subtract (e.g. emptyroom) baseline power spectrum
if isfield(opt, 'baseline')
  powspctrm_baseline = opt.baseline;
  opt = rmfield(opt, 'baseline'); % don't know downstream effects so better remove it
else
  powspctrm_baseline = [];
end
hasbaseline = ~isempty(powspctrm_baseline) && all(size(powspctrm)==size(powspctrm_baseline));

opt.freq_range          = ft_getopt(opt, 'freq_range',        [freqs(find(freqs>0,1,'first')) freqs(end)]);
opt.peak_width_limits   = ft_getopt(opt, 'peak_width_limits', opts_bst.peakwidth.Value{1});
opt.max_peaks           = ft_getopt(opt, 'max_peaks',         opts_bst.maxpeaks.Value{1});
opt.min_peak_height     = ft_getopt(opt, 'min_peak_height',   opts_bst.minpeakheight.Value{1}/10); % convert from dB to B
opt.aperiodic_mode      = ft_getopt(opt, 'aperiodic_mode',    opts_bst.apermode.Value);
opt.peak_threshold      = ft_getopt(opt, 'peak_threshold',    2);   % 2 std dev: parameter for interface simplification
opt.return_spectrum     = ft_getopt(opt, 'return_spectrum',   1);   % SPM/FT: set to 1
opt.border_threshold    = ft_getopt(opt, 'border_threshold',  1);   % 1 std dev: proximity to edge of spectrum, static in Python

% Matlab-only options
opt.power_line          = ft_getopt(opt, 'power_line',        '50'); % for some reason it should be a string, if you don't want a notch, use 'inf'. Brainstorm's default is '60'
opt.peak_type           = ft_getopt(opt, 'peak_type',         opts_bst.peaktype.Value);
opt.proximity_threshold = ft_getopt(opt, 'proximity_threshold', opts_bst.proxthresh.Value{1});
opt.guess_weight        = ft_getopt(opt, 'guess_weight',      opts_bst.guessweight.Value);
opt.thresh_after        = ft_getopt(opt, 'thresh_after',      true);   % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)    

% Output options
opt.sort_type  = opts_bst.sorttype.Value;
opt.sort_param = opts_bst.sortparam.Value;
opt.sort_bands = opts_bst.sortbands.Value;

% Check input frequency bounds
if (any(opt.freq_range <= 0) || opt.freq_range(1) >= opt.freq_range(2))
  ft_error('invalid frequency range for fooof');
end

if exist('fmincon', 'file')
  disp('Using constrained optimization, Guess Weight ignored.')
end

if hasbaseline
  ft_warning('subtracting a baseline powerspctrm estimate before fooof, this is experimental');
  powspctrm = max(powspctrm - powspctrm_baseline, 100*eps(min(powspctrm(:))));
end


% permute is needed because process_fooof expects a singleton second dimension of a 3D matrix...
for k = 1:size(powspctrm,1)
  ft_progress(k./size(powspctrm,1), 'parameterizing spectrum, using FOOOF for trial %d\n', k);
  [fs, fg(k,:)] = process_fooof('FOOOF_matlab', permute(powspctrm(k,:,:), [2 1 3]), freqs, opt, exist('fmincon', 'file'));
end
