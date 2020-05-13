function y = maxstar(x, w)
% maxstar   Log of a sum of exponentials.
%   For vectors, maxstar(x) is equivalent to log(sum(exp(x))).
%   For matrices, maxstar(x) is a row vector and maxstar operates on 
%   each column of x. For N-D arrays, maxstar(x) operates along the
%   first non-singleton dimension.
%
%   maxstar(x,w) is the log of a weighted sum of exponentials,
%   equivalent to log(sum(w.*exp(x))). Vectors w and x must be
%   the same length. For matrix x, the weights w can be input as
%   a matrix the same size as x, or as a vector of the same length 
%   as columns of x. Weights may be zero or negative, but the result
%   sum(w.*exp(x)) must be greater than zero. 
%
%   Acts on first dimensions

if nargin<2
    w = [];
else
    if ~isvector(w)
        error('maxstar: w must be a vector')
    end
    if length(w) ~= size(x,1)
        error('maxstar: weight does not match x')
    end
end
%%
w = w(:);
szx = size(x);
if isempty(w)
    % no weight
    m = max(x);
    y = m + log(sum(exp(bsxfun(@minux,x,m))));
else
    % Move the weight into the exponent xw and find
    % m = max(xw) over terms with positive weights
    wpos = w>0;
    xw = bsxfun(@plus, x(wpos,:), log(w(wpos)));
    m = max(xw);
    exwp = exp( bsxfun(@minus, xw, m) );
    wneg = w<0;
    exwn = exp( x(wneg,:) + bsxfun(@minus,log(-w(wneg)), m) );
    y = m + log(sum(exwp,1) - sum(exwn,1));
end
y = reshape(y,[szx(2:end) 1]);




