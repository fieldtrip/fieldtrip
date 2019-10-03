function [X meaneffective sigmaeffective] = TruncatedGaussian(sigma, range, varargin)
% function X = TruncatedGaussian(sigma, range)
%          X = TruncatedGaussian(sigma, range, n)
%
% Purpose: generate a pseudo-random vector X of size n, X are drawn from
% the truncated Gaussian distribution in a RANGE braket; and satisfies
% std(X)=sigma.
% RANGE is of the form [left,right] defining the braket where X belongs.
% For a scalar input RANGE, the braket is [-RANGE,RANGE].
%
% X = TruncatedGaussian(..., 'double') or
% X = TruncatedGaussian(..., 'single') return an array X of of the
%    specified class.
%
% If input SIGMA is negative, X will be forced to have the same "shape" of
% distribution function than the unbounded Gaussian with standard deviation
% -SIGMA: N(0,-SIGMA). It is similar to calling RANDN and throw away values
% ouside RANGE. In this case, the standard deviation of the truncated
% Gaussian will be different than -SIGMA. The *effective* mean and
% the effective standard deviation can be obtained by calling:
%   [X meaneffective sigmaeffective] = TruncatedGaussian(...)
%
% Example:
% 
% sigma=2;
% range=[-3 5]
% 
% [X meaneff sigmaeff] = TruncatedGaussian(sigma, range, [1 1e6]);
% 
% stdX=std(X);
% fprintf('mean(X)=%g, estimated=%g\n',meaneff, mean(X))
% fprintf('sigma=%g, effective=%g, estimated=%g\n', sigma, sigmaeff, stdX)
% hist(X,64)
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Last update: 19/April/2009

% We keep track this variables so as to avoid calling fzero if
% TruncatedGaussian is called succesively with the same sigma and range
persistent PREVSIGMA PREVRANGE PREVSIGMAC

% shape preserving?
shapeflag = (sigma<0);

% Force inputs to be double class
range = double(range);
if isscalar(range)
    % make sure it's positive
    range=abs(range);
    range=[-range range];
else
    range=sort(range); % right order
end
sigma = abs(double(sigma));

n=varargin;

if shapeflag
    % Prevent the same pdf as with the normal distribution N(0,sigma)
    sigmac = sigma;
else
    if diff(range)^2<12*sigma^2 % This imposes a limit of sigma wrt range
        warning('TruncatedGaussian:RangeSigmaIncompatible', ...
                'TruncatedGaussian: range and sigma are incompatible\n');
        sigmac = Inf;
    elseif isequal([sigma range], [PREVSIGMA PREVRANGE]) % See line #80
        sigmac = PREVSIGMAC; % Do not need to call fzero
    else
        % Search for "sigmac" such that the truncated Gaussian having
        % sigmac in the formula of its pdf gives a standard deviation
        % equal to sigma
        [sigmac res flag] = fzero(@scz,sigma,[],...
                                   sigma^2,range(1),range(2));
        sigmac = abs(sigmac); % Force it to be positive
        if flag<0 % Someting is wrong
            error('TruncatedGaussian:fzerofailled', ...
                  'Could not estimate sigmac\n');
        end
        % Backup the solution
        [PREVSIGMA PREVRANGE PREVSIGMAC] = deal(sigma,range,sigmac);
    end
end

% Compute effective standard deviation
meaneffective=meantrunc(range(1), range(2), sigmac);
sigmaeffective=stdtrunc(range(1), range(2), sigmac);

% Inverse of the cdf functions
if isinf(sigmac)
    % Uniform distribution to maximize the standard deviation within the
    % range. It is like a Gaussian with infinity standard deviation
    if any(strcmpi(n,'single'))
        range = single(range);
    end
    cdfinv = @(y) range(1)+y*diff(range);
else
    c = sqrt(2)*sigmac;
    e = erf(range/c);
    if any(strcmpi(n,'single'))
        % cdfinv will be single class
        c = single(c);
        e = single(e);
    end
    cdfinv = @(y) c*erfinv(e(1)+diff(e)*y);
end

% Generate random variable
X = cdfinv(rand(n{:}));
% Clip to prevent some nasty numerical issues with of erfinv function
% when argument gets close to +/-1
X = max(min(X,range(2)),range(1));

return

end % TruncatedGaussian

function m=meantrunc(lower, upper, s)
% Compute the mean of a trunctated gaussian distribution
if isinf(s)
    m = (upper+lower)/2;
else
    a = (lower/sqrt(2))./s;
    b = (upper/sqrt(2))./s;
    corr = sqrt(2/pi)*(-exp(-b.^2)+exp(-a.^2))./(erf(b)-erf(a));
    m = s.*corr;
end
end % vartrunc

function v=vartrunc(lower, upper, s)
% Compute the variance of a trunctated gaussian distribution
if isinf(s)
    v = (upper-lower)^2/12;
else
    a = (lower/sqrt(2))./s;
    b = (upper/sqrt(2))./s;
    if isinf(a)
        ea=0;
    else
        ea = a.*exp(-a.^2);
    end
    if isinf(b)
        eb = 0;
    else
        eb = b.*exp(-b.^2);
    end
    corr = 1 - (2/sqrt(pi))*(eb-ea)./(erf(b)-erf(a));
    v = s.^2.*corr;
end
end % vartrunc

function stdt=stdtrunc(lower, upper, s)
% Standard deviation of a trunctated gaussian distribution
stdt = sqrt(vartrunc(lower, upper, s)-meantrunc(lower, upper, s).^2);
end % stdtrunc

function res=scz(sc, targetsigma2, lower, upper)
% Gateway for fzero, aim the standard deviation to a target value
res = vartrunc(lower, upper, sc) - targetsigma2 - ...
      meantrunc(lower, upper, sc).^2;
end % scz

% End of file TruncatedGaussian.m
