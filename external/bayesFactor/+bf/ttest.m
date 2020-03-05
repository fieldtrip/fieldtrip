function [bf10,pValue,CI,stats] = ttest(X,varargin)
%TTEST One-sample and paired-sample t-test.
% function [bf10,p,CI,stats] = ttest(X,varargin)    - one sample
% function [bf10,p,CI,stats] = ttest(X,M,varargin)   -one sample,non-zero mean
%
% Calculates Bayes Factor for a one-sample or paired T-test.
%
% X = single sample observations  (a column vector)
% Y = paired observations (column vector) or a scalar mean to compare the samples in X to.
%       [Defaults to 0]
%
% Optional Parm/Value pairs:
% alpha - significance level for the frequentist test. [0.05]
% tail - 'both','right', or 'left' for two or one-tailed tests [both]
% scale - Scale of the Cauchy prior on the effect size  [sqrt(2)/2]
% stats - A struct containing .tstat  , .df , .pvalue .tail and .N - This allows one to
%               calculate BF10 directly from the results of a standard T-Test output.
%
% OUTPUT
% bf10 - The Bayes Factor for the hypothesis that the mean is different
%           from zero
% p - p value of the frequentist hypothesis test
% CI    - Confidence interval for the true mean of X
% stats - Structure with .tstat, .df,
%
% BK - Nov 2018

if isnumeric(X)
    if mod(numel(varargin),2)==0
        % Only X specified
        Y = 0;
        parms = varargin;
    else
        % X and Y specified
        if numel(varargin)>1
            parms = varargin(2:end);
        else
            parms = {};
        end
        Y  =varargin{1};
    end
else
    %Neither X nor Y specified (must be a call with 'stats' specified
    parms = cat(2,X,varargin);
    X=[];Y=[];
end
p=inputParser;
p.addParameter('alpha',0.05);
p.addParameter('tail','both',@(x) (ischar(x)&& ismember(upper(x),{'BOTH','RIGHT','LEFT'})));
p.addParameter('scale',sqrt(2)/2);
p.addParameter('stats',[],@isstruct);
p.parse(parms{:});


if isempty(p.Results.stats)
    % Calculate frequentist from the X and Y data
    tail = p.Results.tail;
    [~,pValue,CI,stats] = ttest(X,Y,'alpha',p.Results.alpha,'tail',tail);
    T = stats.tstat;
    df = stats.df;
    N = numel(X);
else
    % User specified outcome of frequentist test (the builtin ttest), calculate BF from T and
    % df.
    T = p.Results.stats.tstat;
    df = p.Results.stats.df;
    pValue = p.Results.stats.p;
    tail  = p.Results.stats.tail;
    N = p.Results.stats.N;
    CI = [NaN NaN];
end

% Use the formula from Rouder 2009
% This is the formula in that paper, but it does not use the
% scale
% numerator = (1+T.^2/(N-1)).^(-N/2);
% fun  = @(g) ( ((1+N.*g).^-0.5) .* (1+T.^2./((1+N.*g).*(N-1))).^(-N/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );
% This is with scale  (Checked against Morey's R package )
r = p.Results.scale;
numerator = (1+T.^2/df).^(-(df+1)/2);
fun  = @(g) ( ((1+N.*g.*r.^2).^-0.5) .* (1+T.^2./((1+N.*g.*r.^2).*df)).^(-(df+1)/2) .* (2*pi).^(-1/2) .* g.^(-3/2).*exp(-1./(2*g))  );

% Integrate over g
bf01 = numerator/integral(fun,0,inf);
% Return BF10
bf10 = 1./bf01;

switch (tail)
    case 'both'
        % Nothing to do
    case {'left','right'}
        % Adjust the BF using hte p-value as an estimate for the posterior
        % (Morey & Wagenmakers, Stats and Prob Letts. 92 (2014):121-124.
        bf10 = 2*(1-pValue)*bf10;
end
end

