function [idi, rt, rn, bb] = idis(pt, pn, y, z, t, varargin)
%IDIS Integrated Discrimination Improvement between two models
% 
%  Description 
%    [IDI,RT,RN,BB] = IDIS(PT,PN,Y,Z,T,OPTIONS) Returns Integrated
%    Discrimination Improvement (IDI) given two vectors of event probabilities
%    at time T, PT for traditional model and PN for new model, and the
%    observed event or censoring times Y and the corresponding censoring
%    indicators Z (0=event, 1=censored). With these inputs, the estimator 
%    attributed to Pencina et al. in the reference is used, with
%    restriction to time T.
%
%    [IDI,RT,RN,BB] = IDIS(PT,PN,OPTIONS) Returns Integrated
%    Discrimination Improvement (IDI) given two vectors of event probabilities
%    at time T, PT for traditional model and PN for new model. Here, the
%    model-based estimator is used (the "new estimator" in the reference).
%
%    Ouputs RT and RN are the R^2 statistics for the two models. BB are
%    Bayesian bootstrap samples of the IDI distribution.
%
%    OPTIONS is optional parameter-value pair
%      rsubstream - number of a random stream to be used for
%                   simulating dirrand variables. This way same
%                   simulation can be obtained for different models. 
%                   See doc RandStream for more information.
%
%  Reference
%    L. E. Chambless, C. P. Cummiskey, and G. Cui (2011). Several
%    methods to assess improvement in risk prediction models:
%    Extension to survival analysis. Statistics in Medicine
%    30(1):22-38.

% Copyright (C) 2012 Ernesto Ulloa, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
if nargin < 3 || ischar(y)
    % model-based estimator
    model_based_estimator = true;
    
    ip.addRequired('pt', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
    ip.addRequired('pn', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
    ip.addParamValue('rsubstream', 0, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
    if nargin > 2
        % if more than 2 input arguments, there must be four as only a
        % single optional parameter is implemented
        if nargin == 4
            ip.parse(pt, pn, y, z);
        else
            error('Invalid number of arguments.');
        end
    else
        ip.parse(pt, pn); 
    end
else
    % Pencina et al. estimator
    model_based_estimator = false;
    
    ip.addRequired('pt', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
    ip.addRequired('pn', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
    ip.addRequired('y', @(x) isreal(x) && all(isfinite(x(:))))
    ip.addRequired('z', @(x) isreal(x) && all(isfinite(x(:))))
    ip.addRequired('t', @(x) isreal(x) && isscalar(x) && ~isnan(x))
    ip.addParamValue('rsubstream', 0, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
    ip.parse(pt, pn, y, z, t, varargin{:})
end
rsubstream=ip.Results.rsubstream;

if nargout < 4
    % without boostrap
    if model_based_estimator
        rt = rsqr(pt);
        rn = rsqr(pn);
    else
        rt = rsqr(pt, y, z, t);
        rn = rsqr(pn, y, z, t);
    end
else
    % with boostrap
    if model_based_estimator
        if rsubstream == 0
            [rt, bbt] = rsqr(pt);
            [rn, bbn] = rsqr(pn);
        else
            [rt, bbt] = rsqr(pt, 'rsubstream', rsubstream);
            [rn, bbn] = rsqr(pn, 'rsubstream', rsubstream);
        end
    else
        if rsubstream == 0
            [rt, bbt] = rsqr(pt, y, z, t);
            [rn, bbn] = rsqr(pn, y, z, t);
        else
            [rt, bbt] = rsqr(pt, y, z, t, 'rsubstream', rsubstream);
            [rn, bbn] = rsqr(pn, y, z, t, 'rsubstream', rsubstream);
        end
    end
    bb = bbn - bbt;
end

idi = rn - rt;

end




