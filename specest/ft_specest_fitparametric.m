function [F, fs, opt] = ft_specest_fitparametric(powspctrm, freqs, varargin)

% FT_SPECEST_FITPARAMETRIC is a helper function for FT_FREQANALYSIS to fit
% a parametric model to a power spectrum, originally designed to implement
% the technique formerly known as FOOOF (and currently as specparam). It
% can also be used as standalone function.
%
% Use as
%   [F, fs, opt] = ft_specest_fitparametric(powspctrm, freqs[, 'opt', opt, 'method', method])
%
% Input arguments:
%   powspctrm = nrpt x nchan x nfreq power spectrum matrix
%   freqs     = 1 x nfreq vector, frequency axis
%   opt       = struct with user defined options, optional
%   method    = string, (default = 'fooof_matlab'), can also be 'fooof_py',
%                 the latter needs matlab to be able to contact a python interpreter that
%                 has access to fooof
%
% Output arguments:
%   F         = structure-array with output of the fitting procedure
%   fs        = frequency axis
%   opt       = struct with the options as applied during fitting

opt    = ft_getopt(varargin, 'opt',    []);
method = ft_getopt(varargin, 'method', 'fooof_matlab');

if ~contains(method, 'fooof')
  ft_error('methods other than fooof are not supported');
end

% check the size of powspctrm and conform it to how process_fooof expects it, i.e. nrpt x nchan x nfreq
siz = size(powspctrm);
if numel(siz)==2
  powspctrm = shiftdim(powspctrm, -1);
  siz = size(powspctrm);
end
sel = find(siz==numel(freqs));
if ~any(sel==3)
  ft_error('input data matrix should have the frequencies defined along the last axis');
end

% Added experimental functionality to subtract (e.g. emptyroom) baseline power spectrum
if isfield(opt, 'baseline')
  powspctrm_baseline = opt.baseline;
  opt = rmfield(opt, 'baseline'); % don't know downstream effects so better remove it
else
  powspctrm_baseline = [];
end
hasbaseline = ~isempty(powspctrm_baseline) && all(size(powspctrm)==size(powspctrm_baseline));

% the options below are taken from the Brainstorm defaults
opt.freq_range          = ft_getopt(opt, 'freq_range',        [freqs(find(freqs>0,1,'first')) freqs(end)]);
opt.peak_width_limits   = ft_getopt(opt, 'peak_width_limits', [0.5 12]);
opt.max_peaks           = ft_getopt(opt, 'max_peaks',         3);
opt.min_peak_height     = ft_getopt(opt, 'min_peak_height',   3/10);    % convert from dB to B
opt.aperiodic_mode      = ft_getopt(opt, 'aperiodic_mode',    'fixed'); % or fixed
opt.peak_threshold      = ft_getopt(opt, 'peak_threshold',    2);       % 2 std dev: parameter for interface simplification
opt.return_spectrum     = ft_getopt(opt, 'return_spectrum',   1);       % SPM/FT: set to 1
opt.border_threshold    = ft_getopt(opt, 'border_threshold',  1);       % 1 std dev: proximity to edge of spectrum, static in Python

% Matlab-only options, i.e. the python version does not have these options 
opt.power_line          = ft_getopt(opt, 'power_line',        '50');       % for some reason it should be a string, if you don't want notches (which results in a discontinuous frequency axis in the ouptut), use 'inf'. Brainstorm's default is '60'
opt.peak_type           = ft_getopt(opt, 'peak_type',         'gaussian'); % or cauchy, or best
opt.proximity_threshold = ft_getopt(opt, 'proximity_threshold', 2);
opt.guess_weight        = ft_getopt(opt, 'guess_weight',      'none'); % or weak, or strong
opt.thresh_after        = ft_getopt(opt, 'thresh_after',      true);   % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)    
opt.verbose             = ft_getopt(opt, 'verbose',           true);

% Output options
opt.sort_type  = ft_getopt(opt, 'sort_type', 'param'); % or 'band'
opt.sort_param = ft_getopt(opt, 'sort_param', 'frequency'); % or amplitude, or std
opt.sort_bands = ft_getopt(opt, 'sort_bands', {'delta' '2, 4'; 'theta' '5, 7'; 'alpha' '8, 12'; 'beta' '15, 20'; 'gamma1' '30, 59'; 'gamma2' '60, 90'});

% Check input frequency bounds
if (any(opt.freq_range <= 0) || opt.freq_range(1) >= opt.freq_range(2))
  ft_error('invalid frequency range for fooof');
end

switch method
  case {'fooof_matlab' 'fooof'}
    % check for brainstorm functions on the path, and add if needed
    ft_hastoolbox('brainstorm', 1);
    if exist('fmincon', 'file')
      disp('Using constrained optimization, Guess Weight ignored.')
    end

  case 'fooof_py'
    try
      py.fooof.FOOOF;
    catch
      ft_error('not sure whether python version of fooof is accessible from this MATLAB session');
    end
    % initialize FOOOF object
    fm = py.fooof.FOOOF(py.list(num2cell(opt.peak_width_limits)), ...
      opt.max_peaks, ...
      opt.min_peak_height, ...
      opt.peak_threshold, ...
      opt.aperiodic_mode, ...
      opt.verbose);
end

if hasbaseline
  ft_warning('subtracting a baseline powerspctrm estimate before fooof, this is experimental');
  powspctrm = max(powspctrm - powspctrm_baseline, 100*eps(min(powspctrm(:))));
end

% permute is needed because process_fooof expects a singleton second dimension of a 3D matrix...
for k = 1:size(powspctrm,1)
  ft_progress(k./size(powspctrm,1), 'parameterizing spectrum, using FOOOF for trial %d\n', k);
  
  switch method
    case {'fooof_matlab' 'fooof'}
      % this version of the algorithm allegedly loops across channels but requires a permutation of the dimensions (not sure why)
      powspctrm_ = permute(powspctrm(k,:,:), [2 1 3]);
      [fs, F(k,:)] = process_fooof('FOOOF_matlab', powspctrm_, freqs, opt, exist('fmincon', 'file'));
    case 'fooof_py'
      powspctrm_ = shiftdim(powspctrm(k,:,:), 1);

      for i = 1:size(powspctrm_, 1)
        % convert inputs to the required type
        pow_       = py.numpy.array(powspctrm_(i,:));
        freqs_     = py.numpy.array(freqs);
        f_range    = py.list(num2cell(opt.freq_range));

        % run FOOOF fit
        fm.fit(freqs_, pow_, f_range)

        % extract outputs
        res = fm.get_results();
        F(k,i).aperiodic_params = double(res.aperiodic_params);
        F(k,i).peak_params      = double(res.peak_params);
        F(k,i).gaussian_params  = double(res.gaussian_params);
        F(k,i).error            = double(res.error);
        F(k,i).r_squared        = double(res.r_squared);
        F(k,i).power_spectrum   = double(fm.power_spectrum);
        F(k,i).fooofed_spectrum = double(fm.fooofed_spectrum_);
        F(k,i).ap_fit           = double(py.getattr(fm, '_ap_fit'));

        fs = double(fm.freqs);
      end
  end
end
