function [R,edges] = binr(R, nt, nb, opt, par)

% BINR Discretize response matrix.
%
%   -----
%   INPUT
%   -----
%   R   - response matrix to be discretized
%   NT  - number of trials per stimulus
%   NB  - number of bins
%   OPT - binning method option
%   PAR - optional parameter
%
%   -------------------
%   THE RESPONSE MATRIX
%   -------------------
%   L-dimensional responses to S distinct stimuli are stored in a response
%   matrix R of size L-by-T-by-S, T being the maximum number of trials
%   available for any of the stimuli. Emtpy trials, i.e., elements of R not
%   corresponding to a recorded response, can take any value.
%
%   -----------------------------------
%   NUMBER OF TRIALS PER STIMULUS ARRAY
%   -----------------------------------
%   This input specifies the number of trials (responses) recorded for each
%   stimulus. It can be either a scalar (for constant number of trials per
%   stimulus) or an array of length S.
%
%   nt must satisfy the following two conditions:
%   - max(nt) = T
%   - length(nt) = S (only if nt is an array)
%
%   --------------
%   NUMBER OF BINS
%   --------------
%   This value specifies the number of bins used for discretizing the
%   response.
%
%   ------------------------------
%   BINNING OPTIONS AND PARAMETERS
%   ------------------------------
%   The binning option specifies the method used to discretize the response
%   matrix. It can be either any of the options in the table below or the
%   name of a user-defined binning function (see "Building and calling
%   custom binning functions" below). Some of the built-in binning methods
%   allow for one or more parameters to be passed to the function as 
%   specified in the table.
%
%   =======================================================================
%   | OPTION      | DESCRIPTION                                           |
%   =======================================================================
%   | 'eqpop'     | EQUIPOPULATED BINNING                                 |
%   |             | ---------------------                                 |
%   |             | The width of the bins is selected so that each bin    |
%   |             | contains (approximately) the same number of values.   |
%   |             |                                                       |
%   |             | Note: using this option for responses which are not   |
%   |             | continuous in nature (i.e., which contain several     |
%   |             | repeated values) may result in the bins being poorly  |
%   |             | equipopulated.                                        |
%   |-------------|-------------------------------------------------------|
%   | 'eqspace'   | EQUISPACED BINNING                                    |
%   |             | ------------------                                    |
%   |             | The range [a, b], provided as a parameter by the user |
%   |             | in the form of a 2-element array, is divided into bins|
%   |             | of equal width. If no interval is specified, then the |
%   |             | range [m, M] (m and M being the max and min of the    |
%   |             | values to be binned, respectively) is used.           |
%   |             |                                                       |
%   |             | Note: a check is performed on whether the specified   |
%   |             | interval includes all values to be binned.            |
%   |-------------|-------------------------------------------------------|
%   | 'ceqspace'  | CENTERED EQUISPACED BINNING                           |
%   |             | ---------------------------                           |
%   |             | The range [C-D,C+D], C being the mean of the values   |
%   |             | to be binned and D = max(C-m, M-C) (m and M being     |
%   |             | defined as above), is divided into intervals of equal |
%   |             | width.                                                |
%   |-------------|-------------------------------------------------------|
%   | 'gseqspace' | GAUSSIAN EQUISPACED BINNING                           |
%   |             | ---------------------------                           |
%   |             | The range [C-N*STD, C+N*STD] (STD being the standard  |
%   |             | deviation of the values to be binned, C being defined |
%   |             | as above and N being a parameter which is passed to   |
%   |             | the function) is divided into intervals of equal      |
%   |             | width. If no N is specified, N=2 is used.             |
%   |             |                                                       |
%   |             | Note: if any of the values falls outside the selected |
%   |             | range the first and last bin are stretched in order   |
%   |             | to accomodate outliers falling below or above the     |
%   |             | range limits, respectively.                           |
%   =======================================================================
%
%   ---------------------------------------------
%   BUILDING AND CALLING CUSTOM BINNING FUNCTIONS
%   ---------------------------------------------
%
%   Users can define their custom binning methods by means of a binning
%   routine of the form:
%       
%       edges = "func_name"(x, nb)
%
%   where
%
%       "FUNC_NAME" - any valid function name (see Matlab documentation for
%                     instructions on how to buil valid function names)
%       X           - a column array contatining the values which need to
%                     be quantized.
%       NB          - the number of bin that need to be used for
%                     discretization.
%       EDGES       - an (NB+1)-long array of strictly monotonically
%                     increasing values corresponding to the edges of the 
%                     quantization bins
%
%   To call the custom plug-in functons simply pass 'func_name' as a string
%   as the OPT parameter.
%
%   ------
%   OUTPUT
%   ------
%   The fuction returns the response matrix R in which the values have been
%   binned, according to the selected method, into integer ranging from 0
%   to NB.
%
%   See also BINR, ENTROPY, INFORMATION

%   Copyright (C) 2009 Cesare Magri
%   Version 1.0.0

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

if ~ischar(opt)
    msg = 'Third input must be the a string specifying the binning type.';
    error('discretize_R:noEdgesFuncName', msg);
end

switch opt
    case {'eqpop'; 'eqspace'; 'ceqspace'; 'gseqspace';'eqpop_fast';'def_edges'}
        isBuiltinFunc = true;
    otherwise
        isBuiltinFunc = false;
end

edgesFunc = str2func(opt);

[L, maxNt, Ns] = size(R);


if isscalar(nt)
    if nt~=maxNt
        msg = 'max(nt) must be equal to size(R,2).';
        error('binr:nt_vs_size2R', msg);
    else
        mask = ':';
    end
elseif isvector(nt)
    if max(nt)~=maxNt
        msg = 'max(nt) must be equal to size(R,2).';
        error('binr:maxNt_vs_size2R', msg);
    elseif length(nt)~=Ns
        msg = 'size(R,3) must match length(opts.nt). Try transposing nt.';
        error('binr:lengthNt_vs_size3R', msg);
    else
        mask = build_logical_mask(nt, maxNt, Ns);
    end
else
    error('Number of trials per stimulus can be specified either with a scalar or a 1-D array.');
end;

if nargin==5, npar = numel(par); end
edges = nan+zeros(nb+1,L);
for k=1:L
    
    if isBuiltinFunc && nargin==5 && npar==L
        % Built-in function are allowed a more elaborate input:
        edges(:,k) = edgesFunc(R(k,mask(:)), nb, par(k));
    elseif isBuiltinFunc && nargin==5 && npar==L*(nb+1)
        % Built-in function are allowed a more elaborate input:
        edges(:,k) = edgesFunc(R(k,mask(:)), nb, par(:,k));

    elseif isBuiltinFunc && nargin==5
      edges(:,k) = edgesFunc(R(k,mask(:)), nb, par);
    else
        edges(:,k) = edgesFunc(R(k,mask(:)), nb);
        
        if size(edges,1)~=nb+1
            msg = 'Incorrect number of edges.';
            error('binr:incorrectNedges', msg);
        end
    end

    % Performing the actual binning:
    [~, R(k,mask(:))] = histc(R(k,mask(:)), edges(:,k));
end

% Bin values are from 0 to NB:
R = R - 1;