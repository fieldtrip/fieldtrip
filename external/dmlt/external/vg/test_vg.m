clear;
close all;

method = 'dual';%'regression';%

nruns=2;            % number of runs. each run draws new sample dataset
n = 1e2;            % input dimension
ns = round(.05*n);  % sparsity level
betax=1;            % inverse noise in the input
betah=1;            % inverse noise response variance

% number of samples in train, validation and test
p = round(2*n/3);
pv = p;
pt = n;

gamma_mf = zeros(1,nruns);
v_mf=zeros(nruns,n);
m_mf=zeros(nruns,n);
n_mf1 = Inf(1,nruns);
error_mf= Inf(1,nruns);
errorv_mf= Inf(1,nruns);
errort_mf= Inf(1,nruns);
par1_mf= Inf(1,nruns);
par2_mf= Inf(1,nruns);
error_opt = Inf(1,nruns);
errorv_opt = Inf(1,nruns);
errort_opt = Inf(1,nruns);
t_mf= Inf(1,nruns);

for iruns=1:nruns,
    
    fprintf('run = %d\n',iruns);

    % generate input data x (train data), xv (validation data)
    % and xt (test data) and outputs y, yv and yt.
    % input and output data are zero mean.
    % number of samples are p, pv and pt, respectively
    % input dimension is n
    dataok = false;
    while ~dataok

        snonzero=randperm(n);
        snonzero(ns+1:end) = [];
    
        w=sparse(1,n);
        w(snonzero)=1;
        
        % noise response
        sigma=sqrt(1/betah);

        % noise input
        sigmax=sqrt(1/betax);

        % training set
        x=sigmax*randn(n,p);
        x=x-mean(x,2)*ones(1,p);
        dx=sqrt(1/p*sum(x.^2,2));
        x=x./(dx*ones(1,p));
        y=w*x+sigma*randn(1,p);
        y=y-mean(y);
 
        % val set
        xv=sigmax*randn(n,pv);
        xv=xv-mean(xv,2)*ones(1,pv);
        yv=w*xv+sigma*randn(1,pv);
        yv=yv-mean(yv);

        % test set
        xt=sigmax*randn(n,pt);
        xt=xt-mean(xt,2)*ones(1,pt);
        yt=w*xt+sigma*randn(1,pt);
        yt=yt-mean(yt);

        dataok = ~any(isnan(x(:))) && ~any(isnan(xv(:))) && ~any(isnan(xt(:)));
    end
    
    % set arguments for Variational Garrote
    args.x = x; args.xv = xv; args.xt = xt;
    args.y = y; args.yv = yv; args.yt = yt;
    args.method = method;
    
    tic;	
    res(iruns) = vg_method(args);
    t_mf(iruns) = toc;
    
	par1_mf(iruns)=sum(abs(res(iruns).v_mf-w));
	par2_mf(iruns)=sum((res(iruns).v_mf-w).^2);
    
    error_opt(iruns)=1/p*sum((y-w*x).^2,2);
    errorv_opt(iruns)=1/pv*sum((yv-w*xv).^2,2);
    errort_opt(iruns)=1/pt*sum((yt-w*xt).^2,2);
end
    
% results are printed on std output.
fprintf('sigma    = %6.4f\n',sigma); 
fprintf('nruns    = %d\n',nruns); 
fprintf('\topt train error = $%6.4f \\pm %6.4f$\n',mean(error_opt),std(error_opt));
fprintf('\topt val error   = $%6.4f \\pm %6.4f$\n',mean(errorv_opt),std(errorv_opt));
fprintf('\topt test error  = $%6.4f \\pm %6.4f$\n',mean(errort_opt),std(errort_opt));

fprintf('VG RESULTS\n');
fprintf('\tMETHOD = %s\n',method);
fprintf('\tgamma       = $%6.4f \\pm %6.4f$\n',mean([res.gamma_mf]),std([res.gamma_mf]));
fprintf('\tbeta        = $%6.4f \\pm %6.4f$\n',mean([res.beta_mf]),std([res.beta_mf]));
fprintf('\tnon-zero    = $%6.4f \\pm %6.4f$\n',mean([res.n_mf1]),std([res.n_mf1]));
fprintf('\ttrain error = $%6.4f \\pm %6.4f$\n',mean([res.error_mf]),std([res.error_mf]));
fprintf('\tval error   = $%6.4f \\pm %6.4f$\n',mean([res.errorv_mf]),std([res.errorv_mf]));
fprintf('\ttest error  = $%6.4f \\pm %6.4f$\n',mean([res.errort_mf]),std([res.errort_mf]));
fprintf('\tpar1 error  = $%6.4f \\pm %6.4f$\n',mean(par1_mf),std(par1_mf));
fprintf('\tpar2 error  = $%6.4f \\pm %6.4f$\n',mean(par2_mf),std(par2_mf));
fprintf('\tcpu time    = $%6.4f \\pm %6.4f$\n',mean(t_mf),std(t_mf));
