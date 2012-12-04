% Tests Variational Garrote using uncorrelated input

clear;
close all;

method ='regression';%'dual';%

nruns=5;            % number of runs. each run draws new sample dataset
n = 100;            % input dimension
ns = round(.05*n);  % sparsity level
betax=1;            % inverse noise in the input
betah=1;            % inverse noise response variance

% number of samples in training and test sets
p = n;
pt = 2*n;

gamma_mf = zeros(1,nruns);
v_mf=zeros(nruns,n);
m_mf=zeros(nruns,n);
n_mf1 = Inf(1,nruns);
errort_opt = Inf(1,nruns);
errort_mf = Inf(1,nruns);
t_mf= Inf(1,nruns);

for iruns=1:nruns,
    
    fprintf('run = %d\n',iruns);

    % generate input data x (train data), xv (validation data)
    % and xt (test data) and outputs y, yv and yt.
    % input and output data are zero mean.
    % number of samples are p, pv and pt, respectively
    % input dimension is n

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

    % set arguments for Variational Garrote
    args.method = method;

    tic;	
    m = dml.garrote;
    m = m.train(x', y');
    res(iruns) = m.model;
    t_mf(iruns) = toc;

    % test set
    xt=sigmax*randn(n,pt);
    xt=xt-mean(xt,2)*ones(1,pt);
    yt=w*xt+sigma*randn(1,pt);
    yt=yt-mean(yt);

    ytvg = m.test(xt);
    errort_mf(iruns)=1/pt*sum((yt-ytvg).^2,2);

	par1_mf(iruns)=sum(abs(m.model.v_mf-w));
	par2_mf(iruns)=sum((m.model.v_mf-w).^2);
    
    errort_opt(iruns)=1/pt*sum((yt-w*xt).^2,2);
end
    
% results are printed on std output.
fprintf('sigma    = %6.4f\n',sigma); 
fprintf('nruns    = %d\n',nruns); 
fprintf('\topt test error  = $%6.4f \\pm %6.4f$\n',mean(errort_opt),std(errort_opt));

fprintf('VG RESULTS\n');
fprintf('\tMETHOD = %s\n',method);
fprintf('\tgamma       = $%6.4f \\pm %6.4f$\n',mean([res.gamma_mf]),std([res.gamma_mf]));
fprintf('\tbeta        = $%6.4f \\pm %6.4f$\n',mean([res.beta_mf]),std([res.beta_mf]));
fprintf('\tnon-zero    = $%6.4f \\pm %6.4f$\n',mean([res.n_mf1]),std([res.n_mf1]));
fprintf('\ttrain error = $%6.4f \\pm %6.4f$\n',mean([res.error_mf]),std([res.error_mf]));
fprintf('\tval error   = $%6.4f \\pm %6.4f$\n',mean([res.errorv_mf]),std([res.errorv_mf]));
fprintf('\ttest error  = $%6.4f \\pm %6.4f$\n',mean(errort_mf),std(errort_mf));
fprintf('\tpar1 error  = $%6.4f \\pm %6.4f$\n',mean(par1_mf),std(par1_mf));
fprintf('\tpar2 error  = $%6.4f \\pm %6.4f$\n',mean(par2_mf),std(par2_mf));
fprintf('\tcpu time    = $%6.4f \\pm %6.4f$\n',mean(t_mf),std(t_mf));
