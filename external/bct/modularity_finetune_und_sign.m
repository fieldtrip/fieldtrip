function [M q]=modularity_finetune_und_sign(W,qtype,M)
%MODULARITY_FINETUNE_UND_SIGN     Optimal community structure and modularity
%
%   Ci     = modularity_finetune_und_sign(W,[],Ci0);
%   Ci     = modularity_finetune_und_sign(W,'sta',Ci0);
%   [Ci Q] = modularity_finetune_und_sign(W,'sta',Ci0);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges. 
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups. 
%
%   This algorithm is inspired by the Kernighan-Lin fine-tuning algorithm
%   and is designed to refine a previously detected community structure.
%
%   Input:      W,      undirected (weighted or binary) connection matrix
%                       with positive and negative weights
%
%               qtype,  modularity type (see Rubinov and Sporns, 2011)
%                           'sta',  Q_* (default if qtype is not specified)
%                           'pos',  Q_+
%                           'smp',  Q_simple
%                           'gja',  Q_GJA
%                           'neg',  Q_-
%
%               Ci0,    initial community affiliation vector (optional)
%
%
%   Output:     Ci,     refined community affiliation vector
%               Q,      modularity (qtype dependent)
%
%   Note: Ci and Q may vary from run to run, due to heuristics in the
%   algorithm. Consequently, it may be worth to compare multiple runs.
%
%   References:
%   Sun et al. (2008)  Europhysics Lett 86, 28004.
%   Rubinov and Sporns (2011) NeuroImage.
%
%
%   Mika Rubinov, UNSW, 2011

%   Modification History:
%   Mar 2011: Original


n=length(W);                                                %number of nodes/modules
if ~exist('qtype','var') || isempty(qtype);
    qtype = 'sta';
end
if ~exist('M','var') || isempty(M);
    M = 1:n;
else
    [dum dum M] = unique(M(:).');                           %align module indices
end

W0= W.*(W>0);                                               %positive weights matrix
W1=-W.*(W<0);                                               %negative weights matrix
s0=sum(W0(:));                                              %positive sum of weights
s1=sum(W1(:));                                              %negative sum of weights
Knm0=zeros(n,n);                                            %positive node-to-module degree
Knm1=zeros(n,n);                                            %negative node-to-module degree
for m=1:max(M)                                              %loop over modules
    Knm0(:,m)=sum(W0(:,M==m),2);
    Knm1(:,m)=sum(W1(:,M==m),2);
end
Kn0=sum(Knm0,2);                                            %positive node degree
Kn1=sum(Knm1,2);                                            %negative node degree
Km0=sum(Knm0,1);                                            %positive module degree
Km1=sum(Knm1,1);                                            %negative module degree

switch qtype
    case 'smp';  d0 = 1/s0;       d1 = 1/s1;                %dQ = dQ0/s0 - dQ1/s1;
    case 'gja';  d0 = 1/(s0+s1);  d1 = 1/(s0+s1);           %dQ = (dQ0 - dQ1)/(s0+s1);
    case 'sta';  d0 = 1/s0;       d1 = 1/(s0+s1);           %dQ = dQ0/s0 - dQ1/(s0+s1);
    case 'pos';  d0 = 1/s0;       d1 = 0;                   %dQ = dQ0/s0;
    case 'neg';  d0 = 0;          d1 = 1/s1;                %dQ = -dQ1/s1;
    otherwise; error('qtype unknown');
end
if ~s0                                                      %adjust for absent positive weights
    s0=1;
    d0=0;
end
if ~s1                                                      %adjust for absent negative weights
    s1=1;
    d1=0;
end

f=1;                                                        %flag for within-hierarchy search
while f, f=0;
    for u=randperm(n);                                      %loop over all nodes in random order
        ma = M(u);                                          %current module of u
        dQ0 = (Knm0(u,:)+W0(u,u)-Knm0(u,ma)) - Kn0(u).*(Km0+Kn0(u)-Km0(ma))/s0;     %positive dQ
        dQ1 = (Knm1(u,:)+W1(u,u)-Knm1(u,ma)) - Kn1(u).*(Km1+Kn1(u)-Km1(ma))/s1;     %negative dQ
        dQ = d0*dQ0 - d1*dQ1;                               %rescaled changes in modularity
        dQ(ma) = 0;                                         %no changes for same module

        [max_dQ mb] = max(dQ);                              %maximal increase in modularity and corresponding module
        if max_dQ>1e-10, f=1;                               %if maximal increase is positive (equiv. dQ(mb)>dQ(ma))
            M(u) = mb;                                      %reassign module

            Knm0(:,mb)=Knm0(:,mb)+W0(:,u);
            Knm1(:,mb)=Knm1(:,mb)+W1(:,u);
            Knm0(:,ma)=Knm0(:,ma)-W0(:,u);
            Knm1(:,ma)=Knm1(:,ma)-W1(:,u);
            Km0(mb)=Km0(mb)+Kn0(u);
            Km1(mb)=Km1(mb)+Kn1(u);
            Km0(ma)=Km0(ma)-Kn0(u);
            Km1(ma)=Km1(ma)-Kn1(u);
        end
    end
end

[dum dum M]=unique(M(:).');                                 %realign module indices
if nargout==2                                               %compute modularity
    m = M(ones(1,n),:);
    Q0 = (W0-(Kn0*Kn0.')/s0).*(m==m.');
    Q1 = (W1-(Kn1*Kn1.')/s1).*(m==m.');
    q = d0*sum(Q0(:)) - d1*sum(Q1(:));
end