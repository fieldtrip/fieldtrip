sigmoid1 = @(x) 1./(1+exp(-x));

iter=0;
dm=1;
eta=eta0;
eta_w=eta_w0;
m=m+1e-13*randn(1,n);
m=max(m,1e-10);
m=min(m,1-1e-10);

if (beta>=beta_max), beta=beta_max-1; end;
kl_tmp = zeros(1,maxiter);
while ((dm>dmmin)&&(iter<maxiter)&&(beta<beta_max)) ;
	iter=iter+1;
	A=eye(p)+1/p*x'*spdiags(m'./(1-m)'./chi_ii,0,n,n)*x;
    klad=A\y';
	beta=p/(y*klad);
	lambda=beta*klad;
   	mold=m;
	klad1=x*lambda;
   	m=(1-eta)*m+eta*sigmoid1(gamma+(klad1'.^2)./chi_ii'./((1-m).^2)/2/beta/p);
	m=max(m,1e-10);
	m=min(m,1-1e-10);
   	dm=max(abs(m-mold));
	
    kl_tmp(iter)=-p/2*log(beta)-1/2/beta*lambda'*A*lambda+y*lambda - gamma*sum(m)+sum(m.*log(m)+(1-m).*log(1-m));
    if iter>1 
        dkl =(kl_tmp(iter)-kl_tmp(iter-1));
        if dkl>1e-10,
            eta=eta/2;
            if eta < 1e10
                break;
            end
            fprintf('\t[%d] dkl=%e eta=%.3e, dm=%.3e\n',iter,kl_tmp(iter)-kl_tmp(iter-1),eta,dm);
            m=mold;
        end
    end
    
end;
if iter==maxiter, gamma, iter, 
	figure(3)
	plot(1:iter,kl_tmp(1:iter));
	sum(m)
	max(eig(A))
	keyboard;
end;
kl1=-p/2*log(beta)-1/2/beta*lambda'*A*lambda+y*lambda - gamma*sum(m)+sum(m.*log(m)+(1-m).*log(1-m));
v=1./chi_ii'.*m./(1-m).*klad1'/beta/p;

