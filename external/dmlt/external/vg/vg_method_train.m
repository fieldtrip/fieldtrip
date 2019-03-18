function res = vg_method_train(X, Y, args)

% Variational Garrote: performs linear regression with L0-norm penalty (Spike and Slab model)
%   See file test_vg.m for an example of use
%
% required parameters (n input dimension, p samples)
%   X     : n x p (training set, input)
%   Y     : 1 x s (training set, output)
%
% args: optional parameters (default)
%   method      : method for optimization 'dual' or 'regression' for fixed gamma ('dual')
%   maxiter     : maximum number of iterations for optimization for fixed gamma (1e4)
%   max_sum_m   : increases gamma values until sum(m)=max_sum_m  (n/2)
%   beta_max    : increases gamma values until beta=beta_max (1e3)
%   n_gamma     : number of gamma values to scan (50)
%   dmmin       : convergence threshold for mean field error (1e-12)
%   valset      : part of the training set used for validation (0.1*p)

%----------------
% REQUIRED PARAMS

n = size(X,1);

%----------------
% OPTIONAL PARAMS

% method for optimization {dual or regression} for fixed gamma
if ~isfield(args, 'method') method='dual'; else method=args.method; end

% maximum number of iterations for optimization for fixed gamma
if ~isfield(args, 'maxiter') maxiter=1e4; else maxiter=args.maxiter; end

% increases gamma values until sum(m)=max_sum_m
if ~isfield(args, 'max_sum_m') max_sum_m=n/2; else max_sum_m=args.max_sum_m; end

% increases gamma values until beta=beta_max
if ~isfield(args, 'beta_max') beta_max=1e3; else beta_max=args.beta_max; end

% number of gamma values to scan
if ~isfield(args, 'n_gamma') n_gamma=50; else n_gamma=args.n_gamma; end

% convergence threshold for mean field error
if ~isfield(args, 'dmmin') dmmin=1e-12; else dmmin=args.dmmin; end

% part of the training set used for validation (default 0.1*p)
if ~isfield(args, 'valset') valset=ceil(.25*n); else valset=args.valset; end

% randomly split training and validation datasets
dataok = false;
nits = 1;
p = size(X,2)-valset;
pv = valset;
while ~dataok && nits < 10
    it = randperm(size(X,2));
    xv = X(:,it(1:valset));
    yv = Y(:,it(1:valset));
    x = X(:,it(valset+1:end));
    y = Y(:,it(valset+1:end));
    
    % normalize training
    x=x-mean(x,2)*ones(1,p);
    dx=sqrt(1/p*sum(x.^2,2));
    x=x./(dx*ones(1,p));
    y=y-mean(y);
    
    % normalize validation
    xv=xv-mean(xv,2)*ones(1,pv);
    dxv=sqrt(1/pv*sum(xv.^2,2));
    xv=xv./(dxv*ones(1,pv));
    yv=yv-mean(yv);
    
    dataok = ~any(isnan(x(:))) && ~any(isnan(xv(:)));
    nits=nits+1;
end
if nits==10
    error('VG Error: Increase training set size');
end

%----------------
% compute garrote solution for range of gammas
% first from gamma_min to gamma_max and then in
% a second pass from gamma_max to gamma_min.

% C is input data covariance matrix.
if strcmp(method, 'regression')
    if n<=1500,	
        C=x*x'/p;
    end;
end

% b is input output covariance
b=x*y'/p;

% sigma is output variance
sigmay=y*y'/p;

% set gamma range (min, max and step size)
delta=1e-8;
[b2sort,isort]=sort(b.^2,'descend');
bsort=b(isort);
gamma_min=log(delta*sigmay/p/max(abs(b)));
eps_gamma=0.001;
gamma_max=eps_gamma*gamma_min;
gamma_all =linspace(gamma_min,gamma_max,n_gamma);

% initial step size of mean field update
eta0=1; %e-2;
% initial step size for change in w in dual.m
eta_w0=0.02;

% input data variance
chi_ii=1/p*sum(x.^2,2);
if sum(abs(chi_ii-1)>1e-10),
    fprintf('input design matrix is not normalized\n');
    pause
end;

lg=n_gamma;
kl_all=inf(lg,2);
v_all=inf(lg,2,n);
m_all=inf(lg,2,n);
beta_all=inf(lg,2);
v_mf_all=inf(lg,n);
m_mf_all=inf(lg,n);
iter_all=inf(lg,2);
beta_mf_all=inf(1,lg);
error_mf_all=inf(1,lg);
errorv_mf_all=inf(1,lg);
m=zeros(1,n);

% the estimated inverse noise variance beta is initialized as the
% output variance
beta=1/sigmay;
i=0;

% for gamma is gamma_min to gamma_max, or when some criteria are
% satisfied
while (beta<beta_max)&&(i<n_gamma)&&(sum(m)<max_sum_m),
	i=i+1;
	gamma=gamma_all(i);
    eval(method);
	v_all(i,1,:)=v;
	m_all(i,1,:)=m;
	beta_all(i,1)=beta;
	iter_all(i,1)=iter;
	kl_all(i,1)=kl1;
	fprintf('gamma = %f beta = %f sum(m) = %f iter = %d kl = %f\n',gamma,beta,sum(m),iter,kl1);
end;

if beta>=beta_max
 fprintf('-----------------------------------------------------------\n');
 fprintf('beta > beta_max (%.3f > %.3f)\n', beta, beta_max);
 if i>1
    m = squeeze(m_all(i-1,1,:))'; 
 end
end
if sum(m)>=max_sum_m
 fprintf('-----------------------------------------------------------\n');
 fprintf('sum(m) > max_sum_m (%.3f > %.3f)\n', sum(m), max_sum_m);
end

% for gamma is current gamma decreasing to gamma_min 
imax=i-1;
for i=imax:-1:1,
	gamma=gamma_all(i);
    eval(method);
	v_all(i,2,:)=v;
	m_all(i,2,:)=m;
	beta_all(i,2)=beta;
	iter_all(i,2)=iter;
	kl_all(i,2)=kl1;
	fprintf('gamma = %f beta = %f sum(m) = %f iter = %d kl = %f\n',gamma,beta,sum(m),iter,kl1);
end;

% select for each gamma from these two solutions the one with lowest KL 
[klmin, imin]=min(kl_all,[],2);
for i=1:imax,
	v_mf_all(i,:)=squeeze(v_all(i,imin(i),:))';
	m_mf_all(i,:)=squeeze(m_all(i,imin(i),:))';
	beta_mf_all(i)=beta_all(i,imin(i));
	error_mf_all(i)=1/p*sum((y-v_mf_all(i,:)*x).^2,2);
	errorv_mf_all(i)=1/pv*sum((yv-v_mf_all(i,:)*xv).^2,2);
end;

% select the gamma that optimizes the validation error (errorv_mf_all)
[minerrorv i]=min(errorv_mf_all(1:imax));

res.gamma_mf=gamma_all(i);
res.v_mf=v_mf_all(i,:);
res.m_mf=m_mf_all(i,:);
res.n_mf1=sum(res.m_mf>0.5);

res.error_mf=error_mf_all(i);
res.errorv_mf=errorv_mf_all(i);
res.beta_mf=beta_mf_all(i);

if (res.beta_mf==beta_max)
	fprintf('beta_max too small: beta_mf %6.4f, beta_max %6.4f\n', beta_mf(iruns),beta_max);
	pause
end;
if (res.gamma_mf==gamma_min)
	fprintf('gamma at minimum range boundary\n');
end;
if (res.gamma_mf==gamma_max)
	fprintf('gamma at maxium range boundary\n');
end;
