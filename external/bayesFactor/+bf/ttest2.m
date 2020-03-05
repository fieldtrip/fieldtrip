function [bf10,p,CI,stats] = ttest2(X,Y,varargin)
% Calculates Bayes Factor for a two-sample t-test.
% X - Sample 1
% Y - Sample 2 (not necessarily of the same size)
%
% Optional Parm/Value pairs:
% alpha - significance level for the frequentist test. [0.05]
% tail - 'both','right', or 'left' for two or one-tailed tests [both]
%               Note that  'right' means X>Y and 'left' is X<Y
% scale - Scale of the Cauchy prior on the effect size  [sqrt(2)/2]
% stats - A struct containing .tstat  , .df , .pvalue .tail and .N - This allows one to
%               calculate BF10 directly from the results of a standard ttest2 output.
%           Note, however, that .N should be adjusted to nX*nY/(nX+nY) ,
%           and                 .df = nx+ny-2
%           If you call this fuction with data (X, Y) this adjustment is
%           done automatically.
%
% OUTPUT
% bf10 - The Bayes Factor for the hypothesis that the means of the samples are different
% p     - p value of the frequentist hypothesis test
% CI    - Confidence interval for the true mean of X
% stats - Structure with .tstat, .df,resulting from the traditional test.
%
% Internally this code calls bf.ttest for the computation of Bayes Factors.
%
%
% BK - Nov 2018

if isnumeric(X)
    parms = varargin;
else
    %Neither X nor Y specified (this must be a call with 'stats' specified
    parms = cat(2,{X,Y,},varargin);
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
    [~,p,CI,stats] = ttest2(X,Y,'alpha',p.Results.alpha,'tail',tail);
    nX = numel(X);
    nY = numel(Y);
    statsForBf = stats;
    statsForBf.p = p;
    statsForBf.N = nX*nY/(nX+nY);
    statsForBf.df = nX+nY-2;
    statsForBf.tail = tail;
else
    % User specified outcome of frequentist test (the builtin ttest), calculate BF from T and
    % df.
    statsForBf = p.Results.stats;
end

bf10 = bf.ttest('stats',statsForBf);
end
