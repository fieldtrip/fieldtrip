function I = mi_gd_vectorized(x, y, Ym, biascorrect, demeaned, class_means)
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


Ntrl = size(x,1);
Nvar = full(sum(Ym,1));
Nvox = size(Ym,2);
Nclass = size(y,2);
if ~all(Nvar==Nvar(1))
  error('the number of input variables per output sample should  be identical');
end
Nvar = Nvar(1);

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

% allocate memory for class-conditional entropies and covariances
Ntrl_y = sum(y);
Hcond  = zeros(Nvox,Nclass);
Cm     = zeros(Nvox,Nvar,Nvar,Nclass);

% allocate memory for overall entropy and covariance
Hunc = zeros(Nvox,1);
Cx   = zeros(Nvox,Nvar,Nvar);

% input data is class-demeaned, this needs to be accounted for in the
% unconditional entropies
c = diag(sqrt(Ntrl_y))*class_means;%*diag(Ntrl_y);

% ix is needed for bookkeeping which columns go into a single multivariate
% combination
[ix, iy] = find(Ym);

for vi1=1:Nvar
  indx1 = ix(vi1:Nvar:numel(ix));
  x1 = x(:,indx1);
  c1 = c(:,indx1);
  
  Cx(:,vi1,vi1) = sum(x1.^2)+sum(c1.^2);
  for yi=1:Nclass
    tmp = x1(y(:,yi),:);
    Cm(:,vi1,vi1,yi) = sum(tmp.^2);
  end
   
  for vi2=(vi1+1):Nvar
    indx2 = ix(vi2:Nvar:numel(ix));
    x2    = x(:,indx2);
    tmp   = transpose(sum(x1.*x2) + sum(c1.*c(:,indx2)));
    
    Cx(:,vi1,vi2) = tmp;
    Cx(:,vi2,vi1) = tmp;
    for yi=1:Nclass
      tmp   = transpose(sum(x1(y(:,yi),:).*x2(y(:,yi),:)));
      Cm(:,vi1,vi2,yi) = tmp;
      Cm(:,vi2,vi1,yi) = tmp;
    end
  end
end


Cx = Cx / (Ntrl-1);
Cx = cholJM(Cx);
for vi=1:Nvar
  Hunc = Hunc + shiftdim(log(Cx(:,vi,vi)));
end

for yi=1:Nclass
  Cm(:,:,:,yi) = Cm(:,:,:,yi) / (Ntrl_y(yi) - 1);
  Cm(:,:,:,yi) = cholJM(Cm(:,:,:,yi));
  for vi=1:Nvar
    Hcond(:,yi) = Hcond(:,yi)+shiftdim(log(Cm(:,vi,vi,yi)));
  end
end

% apply bias corrections
ln2 = log(2);
if biascorrect

  vars = 1:Nvar;
  
  psiterms_unc = psi((Ntrl - vars)/2) / 2;
  dterm_unc    = (ln2 - log(Ntrl-1)) / 2;
  bias_unc     = Nvar'*dterm_unc + sum(psiterms_unc);
  
  dterm_cond    = (ln2 - log(Ntrl_y-1)) / 2;
  psiterms_cond = zeros(1,Nclass);
  for vi=vars
    idx = (Ntrl_y-vi);
    psiterms_cond = psiterms_cond + psi(idx/2);
  end
  bias_cond = Nvar*dterm_cond + (psiterms_cond/2);
  
  
  Hunc  = Hunc  - bias_unc;
  Hcond = Hcond - ones(Nvox,1)*bias_cond;
end

% class weights
w = Ntrl_y ./ Ntrl;

% compute mutual information
I = Hunc - Hcond*w';

% convert to bits
I = I / ln2;
