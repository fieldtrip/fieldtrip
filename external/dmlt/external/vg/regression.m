sigmoid1 = @(x) 1./(1+exp(-x));

kl_tmp=zeros(1,maxiter);
beta_tmp=zeros(1,maxiter);
dm_tmp=zeros(1,maxiter);
iter=0;
dm=1;
eta=eta0;
m=m+1e-13*randn(1,n);
m=max(m,1e-10);
m=min(m,1-1e-10);
m0 = m;
if (beta>=beta_max), beta=beta_max-1; end;
while ((dm>dmmin)&&(iter<maxiter)&&(beta<beta_max)) 
	iter=iter+1;
	z=(1-m)./m.*diag(C)';
	C1=C+diag(z);
    v=C1\b; 
	w1=v./m';
	betaold=beta;
	beta=1/(sigmay-b'*v);
    
    mold=m;
    m=(1-eta)*m+eta*sigmoid1(gamma+0.5*p*beta*(w1').^2.*diag(C)');
    m=max(m,1e-10);
    m=min(m,1-1e-10);
    dm=max(abs(m-mold));
    kl_tmp(iter)=-p/2*log(beta)+beta*p/2*(v'*C*v+sum((1-m')./m'.*v.^2.*diag(C))-2*b'*v)+beta*p/2*sigmay- gamma*sum(m)+sum(m.*log(m)+(1-m).*log(1-m));
    if iter>1 
        dkl =(kl_tmp(iter)-kl_tmp(iter-1));
        if dkl>1e-10,
            eta=eta/2;
            fprintf('\t[%d] dkl=%e eta=%.3e, dm=%.3e\n',iter,kl_tmp(iter)-kl_tmp(iter-1),eta,dm);
            m=mold;
        end
    end
end;
if iter==maxiter, gamma, iter, 
	figure(3)
	plot(1:iter,kl_tmp(1:iter));
end;
kl1=-p/2*log(beta)+beta*p/2*(v'*C*v+sum((1-m')./m'.*v.^2.*diag(C))-2*b'*v)+beta*p/2*sigmay- gamma*sum(m)+sum(m.*log(m)+(1-m).*log(1-m));
