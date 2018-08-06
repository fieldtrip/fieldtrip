function b = generative_model(A,D,m,modeltype,modelvar,params,epsilon)
%GENERATIVE_MODEL   run generative model code
%
%   B = GENERATIVE_MODEL(A,D,m,modeltype,modelvar,params)
%
%   Generates synthetic networks using the models described in the study by
%   Betzel et al (2016) in Neuroimage.
%
%   Inputs:
%           A,          binary network of seed connections
%           D,          Euclidean distance/fiber length matrix
%           m,          number of connections that should be present in
%                       final synthetic network
%           modeltype,  specifies the generative rule (see below)
%           modelvar,   specifies whether the generative rules are based on
%                       power-law or exponential relationship
%                       ({'powerlaw'}|{'exponential})
%           params,     either a vector (in the case of the geometric
%                       model) or a matrix (for all other models) of
%                       parameters at which the model should be evaluated.
%           epsilon,    the baseline probability of forming a particular
%                       connection (should be a very small number
%                       {default = 1e-5}).
%
%   Output:
%           B,          m x number of networks matrix of connections
%
%
%   Full list of model types:
%   (each model type realizes a different generative rule)
%
%       1.  'sptl'          spatial model
%       2.  'neighbors'     number of common neighbors
%       3.  'matching'      matching index
%       4.  'clu-avg'       average clustering coeff.
%       5.  'clu-min'       minimum clustering coeff.
%       6.  'clu-max'       maximum clustering coeff.
%       7.  'clu-diff'      difference in clustering coeff.
%       8.  'clu-prod'      product of clustering coeff.
%       9.  'deg-avg'       average degree
%       10. 'deg-min'       minimum degree
%       11. 'deg-max'       maximum degree
%       12. 'deg-diff'      difference in degree
%       13. 'deg-prod'      product of degree
%
%
%   Example usage:
%
%       load demo_generative_models_data
%
%       % get number of bi-directional connections
%       m = nnz(A)/2;
% 
%       % get cardinality of network
%       n = length(A);
% 
%       % set model type
%       modeltype = 'neighbors';
% 
%       % set whether the model is based on powerlaw or exponentials
%       modelvar = [{'powerlaw'},{'powerlaw'}];
% 
%       % choose some model parameters
%       params = [-2,0.2; -5,1.2; -1,1.5];
%       nparams = size(params,1);
% 
%       % generate synthetic networks
%       B = generative_model(Aseed,D,m,modeltype,modelvar,params);
%
%       % store them in adjacency matrix format
%       Asynth = zeros(n,n,nparams);
%       for i = 1:nparams; 
%           a = zeros(n); a(B(:,i)) = 1; a = a + a'; 
%           Asynth(:,:,i) = a; 
%       end
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015

if ~exist('epsilon','var')
    epsilon = 1e-5;
end

n = length(D);
nparams = size(params,1);
b = zeros(m,nparams);

switch modeltype
    
    case 'clu-avg'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@plus,clu(:,ones(1,n)),clu')/2;
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_avg(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'clu-diff'
        clu = clustering_coef_bu(A);
        Kseed = abs(bsxfun(@minus,clu(:,ones(1,n)),clu'));
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_diff(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'clu-max'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@max,clu(:,ones(1,n)),clu');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_max(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'clu-min'
        clu = clustering_coef_bu(A);
        Kseed = bsxfun(@min,clu(:,ones(1,n)),clu');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_min(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'clu-prod'
        clu = clustering_coef_bu(A);
        Kseed = clu*clu';
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_clu_prod(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'deg-avg'
        kseed = sum(A,2);
        Kseed = bsxfun(@plus,kseed(:,ones(1,n)),kseed')/2;
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_avg(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'deg-diff'
        kseed = sum(A,2);
        Kseed = abs(bsxfun(@minus,kseed(:,ones(1,n)),kseed'));
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_diff(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'deg-max'
        kseed = sum(A,2);
        Kseed = bsxfun(@max,kseed(:,ones(1,n)),kseed');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_max(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'deg-min'
        kseed = sum(A,2);
        Kseed = bsxfun(@min,kseed(:,ones(1,n)),kseed');
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_min(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'deg-prod'
        kseed = sum(A,2);
        Kseed = (kseed*kseed').*~eye(n);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_deg_prod(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'neighbors'
        Kseed = (A*A).*~eye(n);
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_nghbrs(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'matching'
        Kseed = matching_ind(A);
        Kseed = Kseed + Kseed';
        for iparam = 1:nparams
            eta = params(iparam,1);
            gam = params(iparam,2);
            b(:,iparam) = fcn_matching(A,Kseed,D,m,eta,gam,modelvar,epsilon);
        end
        
    case 'sptl'
        for iparam = 1:nparams
            eta = params(iparam,1);
            b(:,iparam) = fcn_sptl(A,D,m,eta,modelvar{1});
        end
        
end

function b = fcn_clu_avg(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
    K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

    switch mv2
        case 'powerlaw'
            Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
            Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
        case 'exponential'
            Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
            Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
    end
    Ff = Ff.*~A;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_clu_diff(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    K(:,bth) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)')) + epsilon;
    K(bth,:) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)'))' + epsilon;

    switch mv2
        case 'powerlaw'
            Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
            Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
        case 'exponential'
            Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
            Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
    end
    Ff = Ff.*~A;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_clu_max(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    K(:,bth) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
    K(bth,:) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;

    switch mv2
        case 'powerlaw'
            Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
            Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
        case 'exponential'
            Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
            Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
    end
    Ff = Ff.*~A;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_clu_min(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    K(:,bth) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
    K(bth,:) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;

    switch mv2
        case 'powerlaw'
            Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
            Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
        case 'exponential'
            Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
            Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
    end
    Ff = Ff.*~A;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_clu_prod(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end

c = clustering_coef_bu(A);
k = sum(A,2);

Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);

for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    K(bth,:) = (c(bth,:)*c') + epsilon;
    K(:,bth) = (c*c(bth,:)') + epsilon;
    
    switch mv2
        case 'powerlaw'
            Ff(bth,:) = Fd(bth,:).*((K(bth,:)).^gam);
            Ff(:,bth) = Fd(:,bth).*((K(:,bth)).^gam);
        case 'exponential'
            Ff(bth,:) = Fd(bth,:).*exp((K(bth,:))*gam);
            Ff(:,bth) = Fd(:,bth).*exp((K(:,bth))*gam);
    end
    Ff = Ff.*~A;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_deg_avg(A,K,D,m,eta,gam,modelvar,epsilon)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
D = D(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
K = K + epsilon;
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
P = Fd.*Fk(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    switch mv2
        case 'powerlaw'
            Fk(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam;
            Fk(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam);
            Fk(w,:) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam)';
    end
    P = Fd.*Fk(indx);
    b(i) = r;
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_deg_diff(A,K,D,m,eta,gam,modelvar,epsilon)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
D = D(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
K = K + epsilon;
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
P = Fd.*Fk(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    switch mv2
        case 'powerlaw'
            Fk(:,w) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon).^gam;
            Fk(w,:) = ((abs([k - k(w(1)), k - k(w(2))]) + epsilon).^gam)';
        case 'exponential'
            Fk(:,w) = exp((abs([k - k(w(1)), k - k(w(2))]) + epsilon)*gam);
            Fk(w,:) = exp((abs([k - k(w(1)), k - k(w(2))]) + epsilon)*gam)';
    end
    P = Fd.*Fk(indx);
    b(i) = r;
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_deg_min(A,K,D,m,eta,gam,modelvar,epsilon)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
D = D(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
K = K + epsilon;
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
P = Fd.*Fk(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    switch mv2
        case 'powerlaw'
            Fk(:,w) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon].^gam;
            Fk(w,:) = ([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]*gam);
            Fk(w,:) = exp([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]*gam)';
    end
    P = Fd.*Fk(indx);
    b(i) = r;
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_deg_max(A,K,D,m,eta,gam,modelvar,epsilon)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
D = D(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
K = K + epsilon;
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
P = Fd.*Fk(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    switch mv2
        case 'powerlaw'
            Fk(:,w) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon].^gam;
            Fk(w,:) = ([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]*gam);
            Fk(w,:) = exp([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]*gam)';
    end
    P = Fd.*Fk(indx);
    b(i) = r;
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_deg_prod(A,K,D,m,eta,gam,modelvar,epsilon)
n = length(D);
mseed = nnz(A)/2;
k = sum(A,2);
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
D = D(indx);
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
K = K + epsilon;
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
P = Fd.*Fk(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    switch mv2
        case 'powerlaw'
            Fk(:,w) = ([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon].^gam);
            Fk(w,:) = (([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon].^gam)');
        case 'exponential'
            Fk(:,w) = exp([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]*gam);
            Fk(w,:) = exp([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]*gam)';
    end
    P = Fd.*Fk(indx);
    b(i) = r;
    P(b(1:i)) = 0;
end
b = indx(b);

function b = fcn_nghbrs(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
A = A > 0;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
%         gam = abs(gam);
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    x = A(uu,:);
    y = A(:,vv);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K(uu,y) = K(uu,y) + 1;
    K(y,uu) = K(y,uu) + 1;
    K(vv,x) = K(vv,x) + 1;
    K(x,vv) = K(x,vv) + 1;
    switch mv2
        case 'powerlaw'
            Ff(uu,y) = Fd(uu,y).*(K(uu,y).^gam);
            Ff(y,uu) = Ff(uu,y)';
            Ff(vv,x) = Fd(vv,x).*(K(vv,x).^gam);
            Ff(x,vv) = Ff(vv,x)';
        case 'exponential'
            Ff(uu,y) = Fd(uu,y).*exp(K(uu,y)*gam);
            Ff(y,uu) = Ff(uu,y)';
            Ff(vv,x) = Fd(vv,x).*exp(K(vv,x)*gam);
            Ff(x,vv) = Ff(vv,x)';
    end
    Ff(A) = 0;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_matching(A,K,D,m,eta,gam,modelvar,epsilon)
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);
for ii = (mseed + 1):m
    
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    
    updateuu = find(A*A(:,uu));
    updateuu(updateuu == uu) = [];
    updateuu(updateuu == vv) = [];
    
    updatevv = find(A*A(:,vv));
    updatevv(updatevv == uu) = [];
    updatevv(updatevv == vv) = [];
    
    c1 = [A(:,uu)', A(uu,:)];
    for i = 1:length(updateuu)
        j = updateuu(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(uu) = 0;  use(uu+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(uu,j) = epsilon;
            K(j,uu) = epsilon;
        else
            K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,uu) = K(uu,j);
        end
        
    end
    
    c1 = [A(:,vv)', A(vv,:)];
    for i = 1:length(updatevv)
        j = updatevv(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(vv) = 0;  use(vv+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(vv,j) = epsilon;
            K(j,vv) = epsilon;
        else
            K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,vv) = K(vv,j);
        end
    end
    switch mv2
        case 'powerlaw'
            Fk = K.^gam;
        case 'exponential'
            Fk = exp(gam*K);
    end
    Ff = Fd.*Fk.*~A;
    P = Ff(indx);
end
b = find(triu(A,1));

function b = fcn_sptl(A,D,m,eta,modelvar)
n = length(D);
mseed = nnz(A)/2;
switch modelvar
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Fd(indx).*~A(indx);
b = zeros(m,1);
b(1:mseed) = find(A(indx));
for i = (mseed + 1):m
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(i) = r;
    P = Fd(indx);
    P(b(1:i)) = 0;
end
b = indx(b);
