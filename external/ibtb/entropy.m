function varargout = entropy(R, opts, varargin)

%ENTROPY Computes mutual entropy and entropy-like quantities using
% different methods and different bias correction procedures.
%
%   ------
%   SYNTAX
%   ------
%       [...] = entropy(R, opts, output list ...)
%
%   ---------
%   ARGUMENTS
%   ---------
%   R           - Response matrix.
%   opts        - Options structure.
%   output list - List of strings specifying what to compute.
%
%   -------------------
%   THE RESPONSE MATRIX
%   -------------------
%   L-dimensional responses to S distinct stimuli are stored in a response
%   matrix R of size L-by-T-by-S, T being the maximum number of trials
%   available for any of the stimuli. Emtpy trials, i.e., elements of R not
%   corresponding to a recorded response, can take any value.
%
%   ---------------------
%   THE OPTIONS STRUCTURE
%   ---------------------
%   The options structure can include any the following fields:
%
%   opts.nt
%   -------
%       This field specifies the number of trials (responses) recorded for
%       each stimulus. It can be either a scalar (for constant number of
%       trials per stimulus) or an array of length S.
%
%       nt must satisfy the following two conditions:
%       - max(nt) = T
%       - length(nt) = S (if nt is an array)
%
%   opts.method
%   -----------
%       This field specifies which estimation method to use and can be one
%       of the following strings:
%
%       =============================
%       | OPTION  | DECRIPTION      |
%       =============================
%       | 'dr'    | direct method   |
%       | 'gs'    | gaussian method |
%       =============================
%
%       IMPORTANT!
%       ==========
%       The direct method requires the response values to be discretized
%       into non-negative integer values (this is meant only in a numerical
%       sense, the MATLAB variable still needs to be of type double). See
%       function BINR.M for instruction on how to discretize the responses.
%       Failing to properly discretizing the response will result in Matlab
%       crashing.
%
%   opts.bias
%   ---------
%       This field specifies the bias correction procedure. It can be one
%       of the following strings:
%
%       =====================================
%       | OPTION  | DECRIPTION              |
%       =====================================
%       | 'qe'    | quadratic extrapolation |
%       | 'pt'    | Panzeri & Treves 1996   |
%       | 'gsb'   | gaussian bias           |
%       | 'naive' | naive estimates         |
%       =====================================
%
%   opts.btsp (optional)
%   --------------------
%       This field must be a (non-negative) scalar specifying how many 
%       bootstrap estimates to compute.
%
%       Bootstrap estimates are performed by means of pairing stimuli and
%       responses at random and computing the entropy quantities for these
%       random pairings; each estimate corresponds to a different random
%       pairing configuration.
%
%       Computing bootstrap estimates might prove useful, and has been
%       implemented, only for the quantities H(R|S), H_ind(R), H_ind(R|S),
%       Chi(R) and H_sh(R|S).
%
%       See the examples below for additional information on how to use
%       this option.
%
%       DEFAULT: 0.
%
%   opts.verbose (optional)
%   -----------------------
%       If this field exists and is set to true a summary of the selected
%       options is displayed and additional checks are performed on the
%       input variables. No warnings are displayed unless this options is
%       enabled.
%
%       This feature is useful to check whether ENTROPY is being called
%       correctly. It is therefore highly reccomended for new users or when
%       first running of the program with new input options. However, keep
%       in mind that these checks drammatically increases computation time
%       and are thus not reccommended for computationally intensive
%       session.
%
%       DEFAULT: false.
%   
%   ---------------
%   THE OUTPUT LIST
%   ---------------
%   To specify which IT quantities need to compute, one or more of the
%   following strings has to be specified:
%
%       ========================
%       | OPTION  | DECRIPTION |
%       ========================
%       | 'HR'    | H(R)       |
%       | 'HRS'   | H(R|S)     |
%       | 'HlR'   | H_lin(R)   |
%       | 'HiR'   | H_ind(R)   |
%       | 'HiRS'  | H_ind(R|S) |
%       | 'ChiR'  | Chi(R)     |
%       | 'HshR'  | H_sh(R)    |
%       | 'HshRS' | H_sh(R|S)  |
%       ========================
%
%   Outputs are returned IN THE SAME ORDER as that specified in the output
%   list.
%
%   IMPORTANT: Not all combinations of method, bias and output options are
%   available. For example, bias correction 'pt' can only be used with
%   method 'dr'. If 'pt' is called in conjunction with method 'gs', then
%   ENTROPY will apply no bias correction at all (i.e. it will return naive
%   estimates). However, by default, NO message will be displayed. To be
%   prompted with warnings, the verbose option must be set to true (see
%   above). The allowed combinations of method, bias and output options are
%   summarized in the following tables:
%
%       =============================================
%       | DIRECT METHOD                             |
%       =============================================
%       |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%       |-------------------------------------------|
%       | 'HR'    |    X    |   X   |   X   | naive |
%       |---------|---------|-------|-------|-------|
%       | 'HRS'   |    X    |   X   |   X   | naive |
%       |---------|---------|-------|-------|-------|
%       | 'HlR'   |    X    |   X   |   X   | naive |
%       |---------|---------|-------|-------|-------|
%       | 'HiR'   |    X    |   X   | naive | naive |
%       |---------|---------|-------|-------|-------|
%       | 'HiRS'  |    X    |   X   |   X   | naive |
%       |---------|---------|-------|-------|-------|
%       | 'ChiR'  |    X    |   X   | naive | naive |
%       |---------|---------|-------|-------|-------|
%       | 'HshR'  |    X    |   X   |   X   | naive |
%       |---------|---------|-------|-------|-------|
%       | 'HshRS' |    X    |   X   |   X   | naive |
%       =============================================
%
%       =============================================
%       | GAUSSIAN METHOD                           |
%       =============================================
%       |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%       |-------------------------------------------|
%       | 'HR'    |    X    |   X   | naive |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'HRS'   |    X    |   X   | naive |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'HlR'   |    X    |   X   | naive |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'HiR'   |   NaN   |  NaN  |  NaN  |  NaN  |
%       |---------|---------|-------|-------|-------|
%       | 'HiRS'  |    X    |   X   | naive |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'ChiR'  |   NaN   |  NaN  |  NaN  |  NaN  |
%       |---------|---------|-------|-------|-------|
%       | 'HshR'  |    X    |   X   | naive |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'HshRS' |    X    |   X   | naive |   X   |
%       =============================================
%
%   Legend: X    : combination available
%           naive: 'naive' estimate returned
%           NaN  : NaN returned
%
%   --------
%   EXAMPLES
%   --------
%   In the following examples, we assume R to be a 2-by-10-by-3 matrix
%   i.e., R stores 2-dimensional responses to 3 different stimuli. We also
%   assume that, while 10 trials are available for stimulus 1 and 2, only 7
%   trials have been recorded for stimulus 3.
%
%   - Estimate H(R) using all values in R, the direct method and no bias
%     correction:
%
%       opts.nt = [10 10 7];
%       opts.method = 'dr';
%       opts.bias = 'naive';
%       X = entropy(R, opts, 'HR');
%
%   - Compute H_ind(R) and H(R) using the direct method and quandratic
%     extrapolation bias correction:
%
%       opts.nt = [10 10 7];
%       opts.method = 'dr';
%       opts.bias = 'qe';
%       [X, Y] = entropy(R, opts, 'HiR', 'HR');
%
%     where the estimate of H_ind(R) is stored in the X variable and that
%     of H(R) in Y.
%
%   - Compute direct naive gaussian estimtes of H(R) and H(R|S) together 
%     with 20 bootstrap estimates of H(R|S):
%
%       opts.nt = [10 10 7];
%       opts.method = 'gs';
%       opts.bias = 'naive';
%       opts.btsp = 20;
%       [X, Y] = entropy(R, opts, 'HR', 'HRS');
%
%     Note that, in this case, Y is an array of size 21-by-1: Y(1) gives
%     the estimate for H(R|S) computed using the input matrix R; Y(2:21)
%     are are 20 distinct bootstrap estimates of H(R|S).
%
%   -------
%   REMARKS
%   -------
%   - Field-names in the option structure are case-sensitive
%
%   - Ouput options are NOT case sensitive, i.e., they are case-INsensitive
%
%   - It is more efficient to call ENTROPY with several output options
%     rather than calling the function repeatedly. For example:
%
%         [X, Y] = entropy(R, opts, 'HR', 'HRS');
%
%     performs faster than
%
%         X = entropy(R, opts, 'HR');
%         Y = entropy(R, opts, 'HRS');
%
%   - Some MEX files in the toolbox create static arrays which are used
%     to store computations performed in previous calls to the routines.
%     This memory is freed automatically when Matlab is quitted. However,
%     consider using
%
%         clear mex;
%
%     when needing to free all of Matlab's available memory.
%
%   See also BUILD_R_AND_NT, INFORMATION

%   Copyright (C) 2009 Cesare Magri
%   Version: 1.0.4

% -------
% LICENSE
% -------
% This software is distributed free under the condition that:
%
% 1. it shall not be incorporated in software that is subsequently sold;
%
% 2. the authorship of the software shall be acknowledged and the following
%    article shall be properly cited in any publication that uses results
%    generated by the software:
%
%      Magri C, Whittingstall K, Singh V, Logothetis NK, Panzeri S: A
%      toolbox for the fast information analysis of multiple-site LFP, EEG
%      and spike train recordings. BMC Neuroscience 2009 10(1):81;
%
% 3.  this notice shall remain in place in each source file.

% NOTE: Starting with version 1.0.3 the distinction between H_lin(R|S) and
% H_ind(R|S) is not longer made visible to the user. When invoking HiRS the
% program actually computes H_lin(R|S) because of its more desirable
% properties. However, H_ind(R|S) can be still computed according to its
% P_ind(r|s)-based definition by calling the option 'HiRSdef'. The
% distinction between the two quantities is, however, kept throughout the
% program, therefore, remember: HiRS in the program is called HlRS!

persistent previous_pars previous_opts

whereuvar  = strcmpi(varargin, 'uvar');
if sum(whereuvar)
  Runit = varargin{find(whereuvar)+1};
  sel = find(whereuvar);
  sel = [sel sel+1];
  varargin(sel) = [];
end

%persistent pars
if ~isequal(opts,previous_opts)
  % HANDLING INPUTS =========================================================
  pars = build_parameters_structure_v4(R, opts, varargin{:});
  % -------------------------------------------------------------------------
else
	pars = previous_pars;
end

% COMPUTING ENTROPIES =====================================================
Ns    = pars.Ns;
btsp  = pars.btsp;
totNt = sum(pars.Nt);

HRS   = zeros(Ns, btsp+1);
HlRS  = zeros(Ns, btsp+1);
HshRS = zeros(Ns, btsp+1);
ChiR  = zeros(1, btsp+1);
HiR   = zeros(1, btsp+1);

if pars.biasCorrNum==1
    [HR(1), HRS(:,1), HlR(1), HlRS(:,1), HiR(1), HiRS(:,1), ChiR(1), HshR(1), HshRS(:,1)] = quadratic_extrapolation_v2(R, pars);
else
    [HR(1), HRS(:,1), HlR(1), HlRS(:,1), HiR(1), HiRS(:,1), ChiR(1), HshR(1), HshRS(:,1)] = pars.methodFunc(R, pars);
end

% Bootstrap
if any(pars.btsp)
    maxNt = max(pars.Nt);
    
    % Linear indexing of the elements of R filled with trials as specified
    % by NT (not considering the first dimension):
    one2maxNt = (1:maxNt).';
    filledTrialsIndxes    = zeros(maxNt, Ns);
    filledTrialsIndxes(:) = 1:maxNt*Ns;
    filledTrialsIndxes    = filledTrialsIndxes(one2maxNt(:, ones(Ns,1))<=pars.Nt(:, ones(maxNt,1)).');
    
    % No bootstrap is computed for HR, HlR, HshR:
    pars.doHR   = false;
    pars.doHlR  = false;
    pars.doHshR = false;
    
    for k=2:btsp+1
        
        % Randperm (inlining for speed):
        [ignore, randIndxes] = sort(rand(totNt,1));

        % Randomly assigning trials to stimuli as defined by NT:
        R(:, filledTrialsIndxes(randIndxes)) = R(:, filledTrialsIndxes);

        if pars.biasCorrNum==1
            [ignore, HRS(:,k), ignore, HlRS(:,k), HiR(:,k), ignore, ChiR(:,k), ignore, HshRS(:,k)] = ...
                quadratic_extrapolation_v2(R, pars);
        else
            [ignore, HRS(:,k), ignore, HlRS(:,k), HiR(:,k), ignore, ChiR(:,k), ignore, HshRS(:,k)] = ...
                pars.methodFunc(R, pars);
        end
    end
end % ---------------------------------------------------------------------




% ASSIGNING OUTPUTS =======================================================
Ps = pars.Nt ./ totNt;
varargout = cell(pars.Noutput,1);

varargout(pars.HR)    = {HR};
varargout(pars.HRS)   = {sum(  HRS .* Ps(:, ones(btsp+1,1)), 1)};

% In the case Nc = 1 we have:
%   HlR = HiR = ChiR = HshR = HR
%   HlRS = HiRS = HshRS = HR
if pars.Nc==1
    varargout(pars.HlR)   = {HR};
    varargout(pars.HlRS)  = {sum(  HRS .* Ps(:, ones(btsp+1,1)), 1)};
    varargout(pars.HiRS)  = {sum(  HRS .* Ps(:, ones(btsp+1,1)), 1)};
    varargout(pars.HshR)  = {HR};
    varargout(pars.HshRS) = {sum(  HRS .* Ps(:, ones(btsp+1,1)), 1)};
    
    switch pars.methodNum
        case 1
            varargout(pars.HiR)  = {HR};
            varargout(pars.ChiR) = {HR};
            
        % In the gaussian method case we return NaN for consistency with
        % the NC>1 case:
        case 2
            varargout(pars.HiR)  = {NaN};
            varargout(pars.ChiR) = {NaN};
    end
else 
    varargout(pars.HlR)   = {HlR};
    varargout(pars.HlRS)  = {sum( HlRS .* Ps(:, ones(btsp+1,1)), 1)};
    varargout(pars.HiR)   = {HiR};
    varargout(pars.HiRS)  = {sum( HiRS .* Ps, 1)};
    varargout(pars.ChiR)  = {ChiR};
    varargout(pars.HshR)  = {HshR};
    varargout(pars.HshRS) = {sum(HshRS .* Ps(:, ones(btsp+1,1)), 1)};
end

previous_opts = opts;
previous_pars = pars;