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
%                 or 'specparam_py' the latter two need matlab to be able to 
%                 use a python interpreter that has access to
%                 fooof/specparam
%
% Output arguments:
%   F         = structure-array with output of the fitting procedure
%   fs        = frequency axis, this may be a discontinuous one in case
%                 fooof_matlab is used in combination with a power_line option 
%   opt       = struct with the options as applied during fitting
%   pythonversion = path to python with specparam/fooof available, can be
%                empty if pyenv has already been called in the current matlab session 

opt           = ft_getopt(varargin, 'opt',    []);
method        = ft_getopt(varargin, 'method', 'fooof_matlab');
pythonversion = ft_getopt(varargin, 'pythonversion', []);

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

% if undefined, ensure that 0 is not in the freq range
opt.freq_range          = ft_getopt(opt, 'freq_range',        [freqs(find(freqs>0,1,'first')) freqs(end)]);

% these are shared between methods
opt.peak_width_limits   = ft_getopt(opt, 'peak_width_limits', [0.5 12]);
opt.max_n_peaks         = ft_getopt(opt, 'max_n_peaks',       3);
opt.min_peak_height     = ft_getopt(opt, 'min_peak_height',   3/10);    % convert from dB to B
opt.peak_threshold      = ft_getopt(opt, 'peak_threshold',    2);       % 2 std dev: parameter for interface simplification
opt.aperiodic_mode      = ft_getopt(opt, 'aperiodic_mode',    'fixed'); % or fixed
opt.verbose             = ft_getopt(opt, 'verbose',           true);

% Output options, not sure where these are used, probably in fooof_matlab 
opt.sort_type  = ft_getopt(opt, 'sort_type', 'param'); % or 'band'
opt.sort_param = ft_getopt(opt, 'sort_param', 'frequency'); % or amplitude, or std
opt.sort_bands = ft_getopt(opt, 'sort_bands', {'delta' '2, 4'; 'theta' '5, 7'; 'alpha' '8, 12'; 'beta' '15, 20'; 'gamma1' '30, 59'; 'gamma2' '60, 90'});

% Check input frequency bounds
if (any(opt.freq_range <= 0) || opt.freq_range(1) >= opt.freq_range(2))
  ft_error('invalid frequency range for fooof');
end

% add the method specific options, check whether python or brainstorm code is available (if needed), and create the relevant object
switch method
  case {'fooof_matlab' 'fooof'}
    % method specific options adjusted from the Brainstorm defaults
    opt.max_peaks           = opt.max_n_peaks; opt = rmfield(opt, 'max_n_peaks'); % inconsistent naming scheme
    opt.power_line          = ft_getopt(opt, 'power_line',        '50');       % for some reason it should be a string, if you don't want notches (which results in a discontinuous frequency axis in the ouptut), use 'inf'. Brainstorm's default is '60'
    opt.peak_type           = ft_getopt(opt, 'peak_type',         'gaussian'); % or cauchy, or best for brainstorm
    opt.proximity_threshold = ft_getopt(opt, 'proximity_threshold', 2);
    opt.guess_weight        = ft_getopt(opt, 'guess_weight',      'none'); % or weak, or strong
    opt.thresh_after        = ft_getopt(opt, 'thresh_after',      true);   % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)    
    opt.return_spectrum     = ft_getopt(opt, 'return_spectrum',   1);       % SPM/FT: set to 1
    opt.border_threshold    = ft_getopt(opt, 'border_threshold',  1);       % 1 std dev: proximity to edge of spectrum, static in Python

    % check for brainstorm functions on the path, and add if needed
    ft_hastoolbox('brainstorm', 1);
    if exist('fmincon', 'file')
      disp('Using constrained optimization, Guess Weight ignored.')
    end

  case 'fooof_py'
    % check whether a pyenv call has already been executed
    if ~isequal(pyenv().Status, 'Loaded') && isempty(pythonversion)
      ft_error('python is not Loaded, and no pythonversion has been specified');
    elseif ~isequal(pyenv().Status, 'Loaded')
      pyenv('Version', pythonversion, 'ExecutionMode', 'OutOfProcess');
    end

    try
      py.fooof.FOOOF;
    catch
      ft_error('fooof cannot be imported in this MATLAB session; use ''pyenv'' to point MATLAB to a Python installation that has fooof (pip install specparam), see ''doc pyenv''');
    end

    % initialize FOOOF object
    fobj = py.fooof.FOOOF(pyargs( ...
      'peak_width_limits', py.list(num2cell(double(opt.peak_width_limits))), ...
      'max_n_peaks',       double(opt.max_n_peaks), ...
      'min_peak_height',   double(opt.min_peak_height), ...
      'peak_threshold',    double(opt.peak_threshold), ...
      'aperiodic_mode',    char(opt.aperiodic_mode), ...
      'verbose',           logical(opt.verbose)));

  case 'specparam_py'
    % method specific options
    opt.periodic_mode = ft_getopt(opt, 'periodic_mode', 'gaussian'); %or cauchy/skewed_gaussian
    opt.debug         = ft_getopt(opt, 'debug', false);
    opt.verbose       = ft_getopt(opt, 'verbose', true);

    % these are private options for the SpectralFitAlgorithm and may change
    % in future versions: the current code has been written using release
    % candidate 7 of specparam
    opt.ap_percentile_thresh = ft_getopt(opt, 'ap_percentile_thresh', 0.025);
    %opt.ap_guess -> Claude told me this is buggy in rc7, JM hasn't verified this
    opt.ap_bounds            = ft_getopt(opt, 'ap_bounds'); % 2x2 ('fixed'), or 2x3 ('knee'), with lower and upper bounds for the parameters
    opt.cf_bound             = ft_getopt(opt, 'cf_bound', 1.5);
    opt.bw_std_edge          = ft_getopt(opt, 'bw_std_edge', 1);
    opt.gauss_overlap_thresh = ft_getopt(opt, 'gauss_overlap_thresh', 0.75);

    % further options for the curve fitting
    opt.maxfev               = ft_getopt(opt, 'maxfev', 5000);
    opt.tol                  = ft_getopt(opt, 'maxtol', 1e-5);

    % check whether a pyenv call has already been executed
    if ~isequal(pyenv().Status, 'Loaded') && isempty(pythonversion)
      ft_error('python is not Loaded, and no pythonversion has been specified');
    elseif ~isequal(pyenv().Status, 'Loaded')
      pyenv('Version', pythonversion, 'ExecutionMode', 'OutOfProcess');
    end

    try
      py.importlib.import_module('specparam');
    catch
      if opt.debug, rethrow(e); end
      ft_error('specparam cannot be imported in this MATLAB session; use ''pyenv'' to point MATLAB to a Python installation that has specparam (pip install specparam), see ''doc pyenv''');
    end

    args = {'aperiodic_mode',    char(opt.aperiodic_mode), ...
      'periodic_mode',     char(opt.periodic_mode), ...
      'peak_width_limits', py.list(num2cell(double(opt.peak_width_limits))), ...
      'max_n_peaks',       double(opt.max_n_peaks), ...
      'min_peak_height',   double(opt.min_peak_height), ...
      'peak_threshold',    double(opt.peak_threshold), ...
      'ap_percentile_thresh', double(opt.ap_percentile_thresh), ...
      'cf_bound',          double(opt.cf_bound), ...
      'bw_std_edge',       double(opt.bw_std_edge), ...
      'gauss_overlap_thresh', double(opt.gauss_overlap_thresh), ...
      'debug',             logical(opt.debug), ...
      'verbose',           logical(opt.verbose)};
 
    if ~isempty(opt.ap_bounds)
      ap_bounds = py.list({py.list(num2cell(opt.ap_bounds(1,:))), py.list(num2cell(opt.ap_bounds(2,:)))});
      args = [args {'ap_bounds' ap_bounds}];
    end
    % initialize the specparam object, in contrast to FOOOF, the 'aperiodic_mode' 
    % and 'periodic_mode' are arguments of the model object, and 'verbose' and 'debug' too.
    fobj = py.specparam.SpectralModel(pyargs(args{:}));
      

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
      [fs, F(k,:)] = process_fooof('FOOOF_matlab', powspctrm_, freqs, opt, exist('fobjincon', 'file'));
    
    case 'fooof_py'
      
      powspctrm_ = shiftdim(powspctrm(k,:,:), 1);
      freqs_     = py.numpy.asarray(double(freqs(:)'), pyargs('dtype', 'float64')).ravel();
      f_range    = py.list(num2cell(opt.freq_range));

      for i = 1:size(powspctrm_, 1)
        
        % select a spectrum and convert to numpy-array, this is somewhat verbose but hopefully robust
        pow_ = py.numpy.asarray(double(powspctrm_(i,:)), pyargs('dtype','float64')).ravel();

        % run FOOOF fit
        try
          fobj.fit(freqs_, pow_, f_range)
        catch e
          ft_warning('specparam fit failed for repetition %d, channel %d: %s', k, i, e.message);
          continue;
        end

        % extract outputs
        res                     = fobj.get_results();
        F(k,i).aperiodic_params = numpy2matlabdouble(res.aperiodic_params);
        F(k,i).peak_params      = numpy2matlabdouble(res.peak_params);
        F(k,i).gaussian_params  = numpy2matlabdouble(res.gaussian_params);
        F(k,i).error            = numpy2matlabdouble(res.error);
        F(k,i).r_squared        = numpy2matlabdouble(res.r_squared);
        F(k,i).power_spectrum   = numpy2matlabdouble(fobj.power_spectrum);
        F(k,i).fooofed_spectrum = numpy2matlabdouble(fobj.fooofed_spectrum_);
        F(k,i).ap_fit           = numpy2matlabdouble(py.getattr(fobj, '_ap_fit'));

        fs = double(fobj.freqs);
      end

    case 'specparam_py'
       
      powspctrm_ = shiftdim(powspctrm(k,:,:), 1);
      freqs_     = py.numpy.asarray(double(freqs(:)'), pyargs('dtype', 'float64')).ravel();
      f_range    = py.list(num2cell(double(opt.freq_range)));

      for i = 1:size(powspctrm_, 1)
        
        % select a spectrum and convert to numpy-array, this is somewhat verbose but hopefully robust
        pow_ = py.numpy.asarray(double(powspctrm_(i,:)), pyargs('dtype','float64')).ravel();
        
        try
          % the same model object can be re-used, it clears prior results at each fit
          fobj.fit(freqs_, pow_, f_range);
        catch e
            ft_warning('specparam fit failed for repetition %d, channel %d: %s', k, i, e.message);
          continue;
        end
        
        res = fobj.results;
        ap = res.params.aperiodic;
        pe = res.params.periodic;

        F(k,i).aperiodic_params = numpy2matlabdouble(ap.get_params('fit'));
        F(k,i).peak_params      = numpy2matlabdouble(pe.get_params('converted'));
        F(k,i).gaussian_params  = numpy2matlabdouble(pe.get_params('fit'));
        F(k,i).periodic_labels  = dictkeys2matlabcell(pe.indices);
        F(k,i).metrics          = dict2struct(res.metrics.results);
        F(k,i).error            = getthisfield(F(k,i).metrics, 'error_mae');
        F(k,i).r_squared        = getthisfield(F(k,i).metrics, 'gof_rsquared');
        F(k,i).power_spectrum   = numpy2matlabdouble(fobj.data.power_spectrum);
        F(k,i).fooofed_spectrum = numpy2matlabdouble(res.model.get_component('full',      'log'));
        F(k,i).ap_fit           = numpy2matlabdouble(res.model.get_component('aperiodic', 'log'));
        F(k,i).peak_fit         = numpy2matlabdouble(res.model.get_component('peak',      'log'));
        F(k,i).fit_ok           = true;
      end

  end
end

% helper functions for conversion
function y = numpy2matlabdouble(x)

% convert numpy array (or numpy/python scalar) into a MATLAB double, preserving the shape
% This goes through a flat list of floats, which works in all MATLAB versions that
% support Python, in contrast to double(ndarray) which is version dependent for N-D arrays.

x  = py.numpy.asarray(x, pyargs('dtype', 'float64'));
sz = cellfun(@double, cell(x.shape));          % python (C-order) shape
v  = double(py.array.array('d', x.ravel().tolist()));
if isempty(sz)
  y = v(1);                                    % 0-d: scalar
elseif isscalar(sz)
  y = reshape(v, 1, sz);                       % 1-D: row vector
else
  y = permute(reshape(v, fliplr(sz)), numel(sz):-1:1); % N-D: C-order -> column-major
end

function k = dictkeys2matlabcell(d)

% python dict -> cell-array with the keys, in order of the values (column index)
keys = cellfun(@char, cell(py.list(d.keys())), 'UniformOutput', false);
vals = cellfun(@(x) numpy2matlabdouble(x), cell(py.list(d.values())));
[~, ix] = sort(vals);
k = keys(ix);

function s = dict2struct(d)

% python dict of {str: float} -> MATLAB struct
s    = struct;
keys = cell(py.list(d.keys()));
for i = 1:numel(keys)
  key = char(keys{i});
  s.(matlab.lang.makeValidName(key)) = numpy2matlabdouble(d.get(keys{i}));
end

function v = getthisfield(s, name)

if isfield(s, name), v = s.(name); else, v = nan; end
