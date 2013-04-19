function [r2, bb] = rsqr(p, y, z, t, varargin)
%RSQR R^2 statistic given probabilities at time point T.
%
%  Description
%    [R2, BB] = RSQR(P,Y,Z,T,OPTIONS) Returns R^2 statistic given vector
%    P of estimated probabilities of presenting event before time T and
%    its density estimate using Bayesian bootstrap. Y are observed event or 
%    censoring times and Z are the corresponding censoring indicators
%    (0=event, 1=censored). With these inputs, the implemented estimator is
%    R^2 estimator attributed to Pencina et al. in the reference, with 
%    restriction to time T.
%
%    [R2, BB] = RSQR(P,OPTIONS) Returns R^2 statistic given vector
%    P of estimated probabilities of presenting event before time T and
%    its density estimate using Bayesian bootstrap. Here, the
%    model-based estimator is used (the "new estimator" in the reference).
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

% Copyright (C) 2012 Tomi Peltola, Ernesto Ulloa, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
if nargin < 2 || ischar(y)
    % model-based estimator
    model_based_estimator = true;
    
    ip.addRequired('p', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
    ip.addParamValue('rsubstream', 0, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
    if nargin > 1
        % if more than 1 input arguments, there must be 3 as only a
        % single optional parameter is implemented
        if nargin == 3
            ip.parse(p, y, z);
        else
            error('Invalid number of arguments.');
        end
    else
        ip.parse(p);
    end
else
    % Pencina et al. estimator
    model_based_estimator = false;
    
    ip.addRequired('p', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
    ip.addRequired('y', @(x) isreal(x) && all(isfinite(x(:))))
    ip.addRequired('z', @(x) isreal(x) && all(isfinite(x(:))))
    ip.addRequired('t', @(x) isreal(x) && isscalar(x) && ~isnan(x))
    ip.addParamValue('rsubstream', 0, @(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
    ip.parse(p, y, z, t, varargin{:})
end
rsubstream=ip.Results.rsubstream;

if model_based_estimator
    % model-based estimator
    s=1-p;
    Es=mean(s);
    Vs=mean(s.^2)-Es^2;
    r2=Vs/(Es*(1-Es));
    
    if nargout > 1
        % do also bootstrap
        n_replicates = 1000;
        n = length(p);
        
        % get substream
        if rsubstream > 0
            prevstream=setrandstream(0,'mrg32k3a');
            stream.Substream = rsubstream;
        end
        
        qr=dirrand(n, n_replicates);
        bb = zeros(n_replicates, 1);
        
        for i=1:n_replicates
            Esr=wmean(s,qr(:,i));
            Vsr=wmean(s.^2,qr(:,i))-Esr^2;
            bb(i,1)=Vsr/(Esr*(1-Esr));
        end
        
        % set RandStream back to previous
        if rsubstream > 0
            setrandstream(prevstream);
        end
    end
else
    % non-model-based estimator
    event_inds = (y <= t) & (z == 0);
    nonevent_inds = y > t;
    % those censored before or at t (y <= t & z == 1) do not contribute to
    % the statistic as their status is unknown
    
    n_events = sum(event_inds);
    n_nonevents = sum(nonevent_inds);
    if n_events < 1
        error('No events.');
    end
    if n_nonevents < 1
        error('No non-events.');
    end
    
    if nargout < 2
        % no bootstrap
        is = mean(p(event_inds));
        ip = mean(p(nonevent_inds));
        r2 = is - ip;
    else
        % with bootstrapping
        n_replicates = 1000;
        n = length(p);
        
        % get substream
        if rsubstream > 0
            prevstream=setrandstream(0,'mrg32k3a');
            stream.Substream = rsubstream;
        end
        
        % bootstrap weights (first column is unit weights)
        w = [ones(n, 1)/n dirrand(n, n_replicates)];
        
        % set RandStream back to previous
        if rsubstream > 0
            setrandstream(prevstream);
        end
        
        % weights over n sum to one, but we need means over n_events and n_nonevents
        is = (n / n_events) * sum(bsxfun(@times, w(event_inds, :), p(event_inds)));
        ip = (n / n_nonevents) * sum(bsxfun(@times, w(nonevent_inds, :), p(nonevent_inds)));
        bb = is - ip;
        
        r2 = bb(1);
        bb(1) = [];
    end
end

end
