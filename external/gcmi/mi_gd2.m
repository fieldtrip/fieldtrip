function I = mi_gd2(x, y, Ym, biascorrect, demeaned, class_means)
% MI_GD Mutual information (MI) between a Gaussian and a discrete
%         variable in bits
%   I = mi_gd(x,y,Ym) returns the MI between the (possibly multidimensional)
%   Gaussian variable x and the discrete variable y.
%   Rows of x correspond to samples, columns to dimensions/variables. 
%   (Samples first axis)
%   y should contain integer values in the range [0 Ym-1] (inclusive).
%
%   biascorrect : true / false option (default true) which specifies whether
%   bias correction should be applied to the esimtated MI.
%   demeaned : false / true option (default false) which specifies whether the
%   input data already has zero mean (true if it has been copula-normalized)

persistent previous_nvar previous_y bias_unc bias_cond;


% ensure samples first axis for vectors
if isvector(x)
    x = x(:);
end
if ndims(x)~=2
    error('mi_gd: input arrays should be 2d')
end
if isvector(y)
    y = y(:);
else
%    error('mi_gd: only univariate discrete variable supported');
end

Ntrl = size(x,1);
Nvar = size(x,2);

if isequal(previous_y, y)
  computebias = false;
else
  computebias = true;
end

if ~isequal(previous_nvar, Nvar)
  computebias = true;
end

if size(y,1) ~= Ntrl
    error('mi_gd: number of trials do not match');
end

% default option values
if nargin<4
    biascorrect = true;
end
if nargin<5
    demeaned = false;
end

if ~demeaned
    x = bsxfun(@minus,x,sum(x,1)/Ntrl);
end

% class-conditional entropies 
Ntrl_y = sum(y);
Hcond = zeros(1,Ym);
for yi=1:Ym
    xm = x(y(:,yi),:);
    %Ntrl_y(yi) = size(xm,1);
    %xm = bsxfun(@minus,xm,sum(xm,1)/Ntrl_y(yi));
    Cm = (xm'*xm) / (Ntrl_y(yi) - 1);
    chCm = chol(Cm);
    Hcond(yi) = sum(log(diag(chCm)));% + c*Nvar;
end
% class weights
w = Ntrl_y ./ Ntrl;

% input data is class-demeaned, this needs to be accounted for in the
% unconditional entropies
c = diag(sqrt(Ntrl_y))*class_means;%*diag(Ntrl_y);

% unconditional entropy from unconditional Gaussian fit
Cx = (x'*x + c'*c) / (Ntrl-1);
chC = chol(Cx);
Hunc = sum(log(diag(chC)));% + c*Nvar; % the commented out bit drops out in the subtraction below


% apply bias corrections
ln2 = log(2);
if biascorrect
    
    
    if computebias,
      vars = 1:Nvar;
      
      psiterms_unc = psi((Ntrl - vars)/2) / 2;
      dterm_unc    = (ln2 - log(Ntrl-1)) / 2;
      bias_unc     = Nvar*dterm_unc + sum(psiterms_unc);
      
      dterm_cond    = (ln2 - log(Ntrl_y-1)) / 2;
      psiterms_cond = zeros(1,Ym);
      for vi=vars
        idx = (Ntrl_y-vi);
        psiterms_cond = psiterms_cond + psi(idx/2);
      end
      bias_cond = Nvar*dterm_cond + (psiterms_cond/2);
    end
    
    Hunc  = Hunc  - bias_unc;
    Hcond = Hcond - bias_cond;
end

I = Hunc - w*Hcond(:);% sum(w .* Hcond);
% convert to bits
I = I / ln2;

previous_y = y;
previous_nvar = Nvar;
