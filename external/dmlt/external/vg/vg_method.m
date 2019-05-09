function res = vg_method(args)

% Variational Garrote: performs linear regression with L0-norm penalty (Spike and Slab model)
%   See file test_vg.m for an example of use
%
% required parameters (n input dimension, p samples)
%   x     : n x p (training set, input)
%   y     : 1 x s (training set, output)
%   xv    : n x p (validation set, input)
%   yv    : 1 x s (validation set, output)
%
% optional parameters (default)
%   method      : method for optimization 'dual' or 'regression' for fixed gamma ('dual')
%   maxiter     : maximum number of iterations for optimization for fixed gamma (1e4)
%   max_sum_m   : increases gamma values until sum(m)=max_sum_m  (n/2)
%   beta_max    : increases gamma values until beta=beta_max (1e3)
%   n_gamma     : number of gamma values to scan (50)
%   dmmin       : convergence threshold for mean field error (1e-12)

%----------------
% REQUIRED PARAMS
ok = true;
if ~isfield(args, 'x') disp('train input x not provided'); ok = false; else x=args.x; end
if ~isfield(args, 'y') disp('train output y not provided'); ok = false; else y=args.y; end
if ~isfield(args, 'xv') disp('val set input xv not provided'); ok = false; else xv=args.xv; end
if ~isfield(args, 'yv') disp('val set output yv not provided'); ok = false; else yv=args.yv; end
if isfield(args, 'xt') && isfield(args, 'yt')
    xt=args.xt;
    yt=args.yt;
    pt = size(xt,2);
else
    pt = 0;
end

if ~ok 
    return; 
end

n = size(x,1);
p = size(x,2);
pv = size(xv,2);

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
errort_mf_all=inf(1,lg);
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
	if pt>0,
		errort_mf_all(i)=1/pt*sum((yt-v_mf_all(i,:)*xt).^2,2);
	end;
end;

% select the gamma that optimizes the validation error (errorv_mf_all)
[minerrorv i]=min(errorv_mf_all(1:imax));

res.gamma_mf=gamma_all(i);
res.v_mf=v_mf_all(i,:);
res.m_mf=m_mf_all(i,:);
res.n_mf1=sum(res.m_mf>0.5);

res.error_mf=error_mf_all(i);
res.errorv_mf=errorv_mf_all(i);
res.errort_mf=errort_mf_all(i);
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
