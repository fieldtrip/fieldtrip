function [c,bb] = hcs(riskscore,y,z,t,varargin)
%HCS Compute Harrell's C for survival model at given time
%
%  Description
%    [C, BB] = HCS(RISKSCORE,Y,Z,T,OPTIONS) Given a risk score vector
%    RISKSCORE, an observed time vector Y, a censoring indicator column
%    vector Z (0=event, 1=censored) and time T, returns
%    Harrell's C at time T and its estimated density using
%    Bayesian Bootstrap method. Large value of RISKSCORE should predict
%    early event.
%
%    The implemented estimator is called ExtAUC(t)_Count in the reference.
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
ip.addRequired('riskscore',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) isreal(x) && all(isfinite(x(:))))
ip.addRequired('z', @(x) isreal(x) && all(isfinite(x(:))))
ip.addRequired('t', @(x) isreal(x) && isscalar(x) && ~isnan(x))
ip.addParamValue('rsubstream',0,@(x) isreal(x) && isscalar(x) && isfinite(x) && x>0)
ip.parse(riskscore,y,z,t,varargin{:})
rsubstream=ip.Results.rsubstream;


n=size(riskscore,1);

if nargout < 2
    % no bootstrapping
    numer = 0;
    denom = 0;
    for i = 1:n
        % i should have had event before or at t
        if z(i) == 1 || y(i) > t
            continue;
        end
        for j = 1:n
            % require y(i) < y(j) and check if the pair is concordant
            if y(i) >= y(j), continue, end;
            numer = numer + (riskscore(i) > riskscore(j));
            denom = denom + 1;
        end
    end
    c = numer / denom;
else
    % with bootstrapping
    n_replicates = 1000;
    
    % get substream
    if rsubstream > 0
        prevstream=setrandstream(0,'mrg32k3a');
        stream.Substream = rsubstream;
    end
    
    % bootstrap weights (first column is unit weights)
    w = [ones(n, 1) dirrand(n, n_replicates)];
    
    % set RandStream back to previous
    if rsubstream > 0
        setrandstream(prevstream);
    end
    
    % same algorithm as above, now with weights
    numer = zeros(1, size(w, 2));
    denom = zeros(1, size(w, 2));
    for i = 1:n
        if z(i) == 1 || y(i) > t
            continue;
        end
        for j = 1:n
            if y(i) >= y(j), continue, end;
            w_ = w(i, :) .* w(j, :);
            numer = numer + w_ * (riskscore(i) > riskscore(j));
            denom = denom + w_;
        end
    end
    bb = numer ./ denom;
    c = bb(1);  % first replicate is with unit weights
    bb(1) = []; % others are the bootstrap replicates
end

end

