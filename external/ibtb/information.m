function varargout = information(R, opts, varargin)

%INFORMATION Computes mutual information using different methods and
% different bias correction procedures.
%
%   ------
%   SYNTAX
%   ------
%       [...] = information(R, opts, output list...)
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
%       NT must satisfy the following two conditions:
%       - max(nt) = T
%       - length(nt) = S (if nt is an array)
%
%   opts.method
%   -----------
%       This field specifies which estimation method to use and can be one
%       of the following strings:
%
%       --------------------------
%       | 'dr' | Direct method   |
%       | 'gs' | Gaussian method |
%       --------------------------
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
%       -------------------------------------
%       | 'qe'    | Quadratic EXtrapolation |
%       | 'pt'    | Panzeri & Treves 1996   |
%       | 'gsb'   | Gaussian bias           |
%       | 'naive' | Biased naive estimates  |
%       -------------------------------------
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
%       This feature is useful to check whether INFORMATION is being called
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
%       =================================================================
%       | Option  | Description       | Expression (in terms of ENTROPY |
%       |         |                   | output options)                 |
%       =================================================================
%       | 'I'     | I(R;S)            | I     = HR - HRS                |
%       | 'Ish'   | I(R;S) shuffle    | Ish   = HR - HiRS + HshRS - HRS |
%       |---------------------------------------------------------------|
%       | 'IX'    | I(R)              | IX    = HlR - HR                |
%       |---------------------------------------------------------------|
%       | 'ILIN'  | I_lin(S;R)        | ILIN  = HlR - HiRS              |
%       |---------------------------------------------------------------|
%       | 'SYN'   | Syn               | SYN   = HR - HRS - HlR + HiRS   |
%       | 'SYNsh' | Syn shuffle       | SYNsh = HR + HshRS - HRS - HlR  |
%       |---------------------------------------------------------------|
%       | 'ISS'   | I_sig_sim         | ISS   = HiR - HlR               |
%       |---------------------------------------------------------------|
%       | 'IC'    | I_cor             | IC    = HR - HRS + HiRS - HiR   |
%       | 'ICsh'  | I_cor shuffle     | ICsh  = HR + HshRS - HRS - HiR  |
%       |---------------------------------------------------------------|
%       | 'ICI'   | I_cor_ind         | ICI   = ChiR - HiR              |
%       |---------------------------------------------------------------|
%       | 'ICD'   | I_cor_dep         | ICD   = HR - HRS - ChiR + HiRS  |
%       | 'ICDsh' | I_cor_dep shuffle | ICDsh = HR + HshRS - HRS - ChiR |
%       |---------------------------------------------------------------|
%       | 'ILB1'  | I_LB1             | ILB1  = HR - HiRS               |
%       |---------------------------------------------------------------|
%       | 'ILB2'  | I_LB2             | ILB2  = ChiR - HiRS             |
%       =================================================================
%
%   Outputs are returned IN THE SAME ORDER as that specified in the output
%   list.
%
%   IMPORTANT: Not all combinations of method, bias and output options are
%   possible. For example, bias correction 'pt' can only be used together 
%   with method 'dr'. The allowed combinations of method, bias and output
%   options are summarized in the following tables:
%
%       =============================================
%       | DIRECT METHOD                             |
%       =============================================
%       |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%       |-------------------------------------------|
%       | 'I'     |    X    |   X   |   X   |   -   |
%       | 'Ish'   |    X    |   X   |   X   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'IX'    |    X    |   X   |   X   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'ILIN'  |    X    |   X   |   X   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'SYN'   |    X    |   X   |   X   |   -   |
%       | 'SYNsh' |    X    |   X   |   X   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'ISS'   |    X    |   X   |   -   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'IC'    |    X    |   X   |   -   |   -   |
%       | 'ICsh'  |    X    |   X   |   -   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'ICI'   |    X    |   X   |   -   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'ICD'   |    X    |   X   |   -   |   -   |
%       | 'ICDsh' |    X    |   X   |   -   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'ILB1'  |    X    |   X   |   X   |   -   |
%       |---------|---------|-------|-------|-------|
%       | 'ILB2'  |    X    |   X   |   -   |   -   |
%       =============================================
%
%   Legend: X: combination available
%           -: combination NOT permitted
%
%
%       =============================================
%       | GAUSSIAN METHOD                           |
%       =============================================
%       |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%       |-------------------------------------------|
%       | 'I'     |    X    |  n.r. |   -   |   X   |
%       | 'Ish'   |    X    |  n.r. |   -   |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'IX'    |    X    |  n.r. |   -   |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'ILIN'  |    X    |  n.r. |   -   |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'SYN'   |    X    |  n.r. |   -   |   X   |
%       | 'SYNsh' |    X    |  n.r. |   -   |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'ISS'   |   NaN   |  NaN  |   -   |  NaN  |
%       |---------|---------|-------|-------|-------|
%       | 'IC'    |   NaN   |  NaN  |   -   |  NaN  |
%       | 'ICsh'  |   NaN   |  NaN  |   -   |  NaN  |
%       |---------|---------|-------|-------|-------|
%       | 'ICI'   |   NaN   |  NaN  |   -   |  NaN  |
%       |---------|---------|-------|-------|-------|
%       | 'ICD'   |   NaN   |  NaN  |   -   |  NaN  |
%       | 'ICDsh' |   NaN   |  NaN  |   -   |  NaN  |
%       |---------|---------|-------|-------|-------|
%       | 'ILB1'  |    X    |  n.r. |   -   |   X   |
%       |---------|---------|-------|-------|-------|
%       | 'ILB2'  |   NaN   |  NaN  |   -   |  NaN  |
%       =============================================
%
%   Legend: X   : combination available
%           -   : combination NOT permitted
%           n.r.: combination available but not recommended
%           NaN : NaN returned
%
%   --------
%   EXAMPLES
%   --------
%   In the following examples, we assume R to be a 2-by-10-by-3 matrix
%   i.e., R stores 2-dimensional responses to 3 different stimuli. We also
%   assume that, while 10 trials are available for stimulus 1 and 2, only 7
%   trials have been recorded for stimulus 3.
%
%   - Estimate I(S;R) using the direct method and no bias corrections
%
%       opts.nt = [10 10 7];
%       opts.method = 'dr';
%       opts.bias = 'naive';
%       X = information(R, opts, 'I');
%
%   - Estimate I(S;R) and I_shuffle(S;R) using direct method and the
%     quadratic extrapolation bias correction
%
%       opts.nt = [10 10 7];
%       opts.method = 'dr';
%       opts.bias = 'gsb';
%       [X, Y] = information(R, opts, 'Ish', 'I');
%
%     where the estimate of I_shuffle(S;R) is stored in the X and that of
%     I(S;R) in Y.
%
%   - Compute gaussian naive estimate of I(S;R) together with 20 bootstrap
%     estimates:
%
%       opts.nt = [10 10 7];
%       opts.method = 'gs';
%       opts.bias = 'naive';
%       opts.btsp = 20;
%       X = information(R, opts, 'I');
%
%     Note that, in this case, Y is an array of size 21-by-1: Y(1) gives
%     the estimate for I(S;R) computed using the input matrix R; Y(2:21)
%     are are 20 distinct bootstrap estimates of I(S;R).
%
%   -------
%   REMARKS
%   -------
%   - Field-names in the option structure are case-sensitive
%
%   - Ouput options are case INsensitive
%
%   - It is more efficient to call INFORMATION with several output options
%     rather than calling the function repeatedly. For example:
%
%         [X, Y] = information(R, opts, 'I', 'Ish');
%
%     is faster than
%
%         X = information(R, opts, 'I');
%         Y = information(R, opts, 'Ish');
%
%   - Some MEX files in the toolbox create static arrays which are used
%     to store computations performed in previous calls to the routines.
%     This memory is freed automatically when Matlab is quitted. However,
%     consider using
%
%         clear mex;
%
%     when needing to free all of Matlab's available memory.

% NOTE:
% The function also computes IXS and IXSsh according to the following
% equations:
%
%       =================================================================
%       | Option  | Description       | Expression (in terms of ENTROPY |
%       |         |                   | output options)                 |
%       =================================================================
%       | 'IXS'   | I(R|S)            | IXS   = HiRS - HRS              |
%       | 'IXSsh' | I(R|S) shuffle    | IXSsh = HshRS - HRS             |
%       =================================================================
%
% However, since this quantity is meaningful only for L=2, this feature
% is kept hidden. The combination tables for this quantity are
% as follows:
%
%       =============================================
%       | DIRECT METHOD                             |
%       =============================================
%       |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%       |-------------------------------------------|
%       | 'IXS'   |    X    |   X   |   X   |   -   |
%       | 'IXSsh' |    X    |   X   |   X   |   -   |
%       =============================================
%
%   Legend: X: combination available
%           -: combination NOT available
%
%
%       =============================================
%       | GAUSSIAN METHOD                           |
%       =============================================
%       |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%       |-------------------------------------------|
%       | 'IXS'   |    X    |  n.r. |   -   |   X   |
%       | 'IXSsh' |    X    |  n.r. |   -   |   X   |
%       =============================================
%
%   Legend: X   : combination available
%           -   : combination NOT available
%           n.r.: combination available but not recommended

%   Copyright (C) 2009 Cesare Magri
%   Version: 1.0.5

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

% NOTE:
% In this function HiRS is indeed output option 'HiRS' of function
% entropy.m, thus corresponding to the variable HlRS in the same function.

whereI     = strcmpi(varargin, 'i');
whereIsh   = strcmpi(varargin, 'ish');
whereIX    = strcmpi(varargin, 'ix');
whereIXS   = strcmpi(varargin, 'ixs');
whereIXSsh = strcmpi(varargin, 'ixssh');
whereILIN  = strcmpi(varargin, 'ilin');
whereSYN   = strcmpi(varargin, 'syn');
whereSYNsh = strcmpi(varargin, 'synsh');
whereISS   = strcmpi(varargin, 'iss');
whereIC    = strcmpi(varargin, 'ic');
whereICsh  = strcmpi(varargin, 'icsh');
whereICI   = strcmpi(varargin, 'ici');
whereICD   = strcmpi(varargin, 'icd');
whereICDsh = strcmpi(varargin, 'icdsh');
whereILB1  = strcmpi(varargin, 'ilb1');
whereILB2  = strcmpi(varargin, 'ilb2');

whereuvar  = strcmpi(varargin, 'uvar');
if sum(whereuvar)
  Runit = varargin{find(whereuvar)+1};
  sel = find(whereuvar);
  sel = [sel sel+1];
  varargin(sel) = [];
end

doI     = any(whereI);
doIsh   = any(whereIsh);
doIX    = any(whereIX);
doIXS   = any(whereIXS);
doIXSsh = any(whereIXSsh);
doILIN  = any(whereILIN);
doSYN   = any(whereSYN);
doSYNsh = any(whereSYNsh);
doISS   = any(whereISS);
doIC    = any(whereIC);
doICsh  = any(whereICsh);
doICI   = any(whereICI);
doICD   = any(whereICD);
doICDsh = any(whereICDsh);
doILB1  = any(whereILB1);
doILB2  = any(whereILB2);

% Checks ------------------------------------------------------------------
specifiedOutputOptsVec = ...
    [doI doIsh doIX doIXS doIXSsh doILIN doSYN doSYNsh doISS doIC doICsh doICI doICD doICDsh doILB1 doILB2];
NspecifiedOutputOpts = sum(specifiedOutputOptsVec);
lengthVarargin = length(varargin);
if NspecifiedOutputOpts~=lengthVarargin
    msg = 'Unknown selection or repeated option in output list.';
    error('information:unknownOutputOpt', msg);
end

% Restrictions on possible combinantions ----------------------------------
if strcmpi(opts.method, 'dr')
    % Can't apply bias-correction gsb with method dr:
    if strcmpi(opts.bias, 'gsb')
        msg = 'Bias correction ''gsb'' can only be used in conjunction with method ''gs''.';
        error('Information:drMethodAndGsbBias', msg);
    end

    % Can't compute ISS, IC, ICsh, ICI, ICD, ICDsh or ILB2 for
    % bias-correction pt
    if strcmpi(opts.bias, 'pt') && (doISS || doIC || doICsh || doICI || doICD || doICDsh || doILB2)
        msg = 'One or more of the selected output options are not available for bias correction ''pt''.';
        error('information:ptBiasAndNonAvailableOutputOpt', msg);
    end
end

% Default verbose value:
if ~isfield(opts, 'verbose')
    opts.verbose = false;
end;

isGaussianMethod = false;
if strcmpi(opts.method, 'gs')
    isGaussianMethod = true;
    
    % Gaussian and QE is not recommended:
    if opts.verbose && strcmpi(opts.bias, 'qe')
        msg = 'Usage of bias correction ''qe'' in conjunction with gaussian method is not recommended.';
        warning('Information:gsMethodAndQeBias', msg);
    end
    
    % Gaussian and PT is not allowed:
    if strcmpi(opts.bias, 'pt')
        msg = 'Bias correction ''pt'' can only be used in conjunction with method ''dr''.';
        error('Information:gsMethodAndPtBias', msg);
    end
end

allOutputOpts = {'HR' 'HRS' 'HlR' 'HiR' 'HiRS' 'ChiR' 'HshRS'};
positionInOuputOptsList = 0;

% What needs to be computed by ENTTROPY -----------------------------------
% Need to compute H(R)? 
doHR = false;
if doI || doIsh || doIX || doSYN || doSYNsh || doIC || doICsh || doICD || doICDsh || doILB1
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHR = positionInOuputOptsList;
    doHR = true;
end

% Need to compute H(R|S)?
doHRS = false;
if doI || doIsh || doIXS || doIXSsh || doSYN || doSYNsh || doIC || doICsh || doICD || doICDsh
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHRS = positionInOuputOptsList;
    doHRS = true;
end

% Need to compute H_lin(R)?
doHlR = false;
if doIX || doSYN || doILIN || doSYNsh || doISS
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHlR = positionInOuputOptsList;
    doHlR = true;
end

% Need to compute H_ind(R)?
doHiR = false;
if (doISS || doIC || doICsh || doICI) && ~isGaussianMethod
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHiR = positionInOuputOptsList;
    doHiR = true;
end

% Need to compute H_ind(R|S)?
doHiRS = false;
if doIsh || doIXS || doILIN || doSYN || doIC || doICD || doILB1 || doILB2
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHiRS = positionInOuputOptsList;
    doHiRS = true;
end

% Need to compute Chi(R)?
doChiR = false;
if (doICI || doICD || doICDsh || doILB2) && ~isGaussianMethod
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereChiR = positionInOuputOptsList;
    doChiR = true;
end

% Need to compute H_sh(R|S)?
doHshRS = false;
if doIsh || doIXSsh || doSYNsh || doICsh || doICDsh
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHshRS = positionInOuputOptsList;
    doHshRS = true;
end

% Computing information theoretic quantities ------------------------------
outputOptsList = allOutputOpts([doHR doHRS doHlR doHiR doHiRS doChiR doHshRS]);

H = cell(positionInOuputOptsList, 1);
if ~exist('Runit', 'var')
  [H{:}] = entropy(R, opts, outputOptsList{:});
else
  [H{:}] = entropy(R, opts, outputOptsList{:}, 'uvar', Runit);
end
  
% Assigning output --------------------------------------------------------
varargout = cell(length(varargin),1);

% I = HR - HRS
if doI
    varargout(whereI) = {H{whereHR} - H{whereHRS}};
end

% Ish = HR - HiRS + HshRS - HRS
if doIsh
    varargout(whereIsh) = {H{whereHR} - H{whereHiRS} + H{whereHshRS} - H{whereHRS}};
end

% IX = HlR - HR
if doIX
    varargout(whereIX) = {H{whereHlR} - H{whereHR}};
end

% IXS = HiRS - HRS
if doIXS
    varargout(whereIXS) = {H{whereHiRS} - H{whereHRS}};
end

% IXSsh = HshRS - HRS
if doIXSsh
    varargout(whereIXSsh) = {H{whereHshRS} - H{whereHRS}};
end

% ILIN = HlR - HiRS
if doILIN
    varargout(whereILIN) = {H{whereHlR} - H{whereHiRS}};
end

% SYN = HR - HRS - HlR + HiRS
if doSYN
    varargout(whereSYN) = {H{whereHR} - H{whereHRS} - H{whereHlR} + H{whereHiRS}};
end

% SYNsh = HR + HshRS - HRS - HlR
if doSYNsh
    varargout(whereSYNsh) = {H{whereHR} + H{whereHshRS} - H{whereHRS} - H{whereHlR}};
end

% ISS = HiR - HlR
if doISS
    if ~isGaussianMethod
        varargout(whereISS) = {H{whereHiR} - H{whereHlR}};
    else
        varargout(whereISS) = {NaN};
    end
end

% IC = HR - HRS + HiRS - HiR
if doIC
    if ~isGaussianMethod
        varargout(whereIC) = {H{whereHR} - H{whereHRS} + H{whereHiRS} - H{whereHiR}};
    else
        varargout(whereIC) = {NaN};
    end
end

% ICsh  = HR + HshRS - HRS - HiR
if doICsh
    if ~isGaussianMethod
        varargout(whereICsh) = {H{whereHR} + H{whereHshRS} - H{whereHRS} - H{whereHiR}};
    else
        varargout(whereICsh) = {NaN};
    end
end

% ICI= ChiR - HiR
if doICI
    if ~isGaussianMethod
        varargout(whereICI) = {H{whereChiR} - H{whereHiR}};
    else
        varargout(whereICI) = {NaN};
    end
end

% ICD = HR - HRS - ChiR + HiRS
if doICD
    if ~isGaussianMethod
        varargout(whereICD) = {H{whereHR} - H{whereHRS} - H{whereChiR} + H{whereHiRS}};
    else
        varargout(whereICD) = {NaN};
    end
end

% ICDsh = HR + HshRS - HRS - ChiR
if doICDsh
    if ~isGaussianMethod
        varargout(whereICDsh) = {H{whereHR} + H{whereHshRS} - H{whereHRS} - H{whereChiR}};
    else
        varargout(whereICDsh) = {NaN};
    end
end

% ILB1 = HR - HiRS
if doILB1
    varargout(whereILB1) = {H{whereHR} - H{whereHiRS}};
end

% ILB2 = ChiR - HiRS
if doILB2
    if ~isGaussianMethod
        varargout(whereILB2) = {H{whereChiR} - H{whereHiRS}};
    else
        varargout(whereILB2) = {NaN};
    end
end