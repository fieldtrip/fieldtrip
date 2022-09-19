function H = ent_g(x, biascorrect)
% ENT_G Entropy of a Gaussian variable in bits
%   H = ent_g(x) returns the entropy of a (possibly 
%   multidimensional) Gaussian variable x with bias correction.
%   Rows of x correspond to samples, columns to dimensions/variables. 
%   (Samples first axis)
%
%   biascorrect : true / false option (default true) which specifies
%   whether bias correction should be applied to the estimated entropy.

if isvector(x)
    x = x(:);
end
if ndims(x)~=2
    error('ent_g: input arrays should be 2d')
end
[Ntrl, Nvarx] = size(x);

if nargin < 2
    % default is to apply bias correction
    biascorrect = true;
end

% demean data
gx.m = sum(x,1)/Ntrl;
x = bsxfun(@minus,x,gx.m);

% covariance
C = (x'*x) / (Ntrl - 1);
chC = chol(C);

% entropy in nats
HX = sum(log(diag(chC))) + 0.5*Nvarx*(log(2*pi)+1);

ln2 = log(2);
if biascorrect
    psiterms = psi((Ntrl - (1:Nvarx))/2) / 2;
    dterm = (ln2 - log(Ntrl-1)) / 2;
    HX = (HX - Nvarx*dterm - sum(psiterms));
end

% convert to bits
H = HX / ln2;
