function [Ci Q] = modularity_louvain_und_sign(W,qtype)
%MODULARITY_LOUVAIN_UND_SIGN     Optimal community structure and modularity
%
%   Ci     = modularity_louvain_und_sign(W);
%   Ci     = modularity_louvain_und_sign(W,'sta');
%   [Ci Q] = modularity_louvain_und_sign(W,'sta');
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges. 
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups. 
%
%   The Louvain algorithm is a fast and accurate community detection 
%   algorithm (at the time of writing).
%
%   Input:      W       undirected (weighted or binary) connection matrix
%                       with positive and negative weights
%
%               qtype,  modularity type (see Rubinov and Sporns, 2011)
%                           'sta',  Q_* (default if qtype is not specified)
%                           'pos',  Q_+
%                           'smp',  Q_simple
%                           'gja',  Q_GJA
%                           'neg',  Q_-
%
%   Output:     Ci,     community affiliation vector
%               Q,      modularity (qtype dependent)
%
%   Note: Ci and Q may vary from run to run, due to heuristics in the
%   algorithm. Consequently, it may be worth to compare multiple runs.
%
%   References:
%   Rubinov and Sporns (2011) NeuroImage.
%   Blondel et al. (2008)  J. Stat. Mech. P10008.
%
%
%   Mika Rubinov, UNSW, 2011

%   Modification History:
%   Mar 2011: Original


N=length(W);                                        %number of nodes
if ~exist('qtype','var') || isempty(qtype);
    qtype = 'sta';
end

W0=W.*(W>0);                                        %positive weights matrix
W1=-W.*(W<0);                                       %negative weights matrix
s0=sum(W0(:));                                      %weight of positive links
s1=sum(W1(:));                                      %weight of negative links

switch qtype
    case 'smp';  d0 = 1/s0;       d1 = 1/s1;        %dQ = dQ0/s0 - dQ1/s1;
    case 'gja';  d0 = 1/(s0+s1);  d1 = 1/(s0+s1);   %dQ = (dQ0 - dQ1)/(s0+s1);
    case 'sta';  d0 = 1/s0;       d1 = 1/(s0+s1);   %dQ = dQ0/s0 - dQ1/(s0+s1);
    case 'pos';  d0 = 1/s0;       d1 = 0;           %dQ = dQ0/s0;
    case 'neg';  d0 = 0;          d1 = 1/s1;        %dQ = -dQ1/s1;
    otherwise; error('qtype unknown');
end
if ~s0                                              %adjust for absent positive weights
    s0=1;
    d0=0;
end
if ~s1                                              %adjust for absent negative weights
    s1=1;
    d1=0;
end

h=2;                                                %hierarchy index
n=N;                                                %number of nodes in hierarchy
Ci={[],1:n};                                        %hierarchical module assignments
Q={-1,0};                                           %hierarchical modularity values
while Q{h}-Q{h-1}>1e-10
    Kn0=sum(W0);                                    %positive node degree
    Kn1=sum(W1);                                    %negative node degree
    Km0=Kn0;                                        %positive module degree
    Km1=Kn1;                                        %negative module degree
    Knm0=W0;                                        %positive node-to-module degree
    Knm1=W1;                                        %negative node-to-module degree

    M=1:n;                                          %initial module assignments

    f=1;                                            %flag for within-hierarchy search
    while f, f=0;
        for u=randperm(n);                          %loop over all nodes in random order
            ma = M(u);                              %current module of u
            dQ0 = (Knm0(u,:)+W0(u,u)-Knm0(u,ma)) - Kn0(u).*(Km0+Kn0(u)-Km0(ma))/s0;     %positive dQ
            dQ1 = (Knm1(u,:)+W1(u,u)-Knm1(u,ma)) - Kn1(u).*(Km1+Kn1(u)-Km1(ma))/s1;     %negative dQ
            dQ = d0*dQ0 - d1*dQ1;                   %rescaled changes in modularity
            dQ(ma) = 0;                             %no changes for same module
            
            [max_dQ mb] = max(dQ);                  %maximal increase in modularity and corresponding module
            if max_dQ>1e-10, f=1;                   %if maximal increase is positive (equiv. dQ(mb)>dQ(ma))
                M(u) = mb;                          %reassign module

                Knm0(:,mb)=Knm0(:,mb)+W0(:,u);      %change positive node-to-module degrees
                Knm0(:,ma)=Knm0(:,ma)-W0(:,u);
                Knm1(:,mb)=Knm1(:,mb)+W1(:,u);      %change negative node-to-module degrees
                Knm1(:,ma)=Knm1(:,ma)-W1(:,u);
                Km0(mb)=Km0(mb)+Kn0(u);             %change positive module degrees
                Km0(ma)=Km0(ma)-Kn0(u);
                Km1(mb)=Km1(mb)+Kn1(u);             %change negative module degrees
                Km1(ma)=Km1(ma)-Kn1(u);
            end
        end
    end

    h=h+1;
    Ci{h}=zeros(1,N);
    [dum dum M]=unique(M);                          %realign module indices
    for u=1:n                                       %loop through initial module assignments
        Ci{h}(Ci{h-1}==u)=M(u);                     %assign new modules
    end

    n=max(M);                                       %number of new nodes
    w0=zeros(n);                                    %new positive weights matrix
    w1=zeros(n);                                    %new negative weights matrix
    for u=1:n
        for v=u:n
            w0(u,v)=sum(sum(W0(M==u,M==v)));        %pool positive weights of nodes in same module
            w1(u,v)=sum(sum(W1(M==u,M==v)));        %pool negative weights of nodes in same module
            w0(v,u)=w0(u,v);
            w1(v,u)=w1(u,v);
        end
    end
    W0=w0;
    W1=w1;

    %compute modularity
    Q0 = sum(diag(W0)) - sum(sum(W0^2))/s0;         %contribution of positive weights
    Q1 = sum(diag(W1)) - sum(sum(W1^2))/s1;         %contribution of negative weights
    Q{h} = d0*Q0 - d1*Q1;
end

Ci=Ci{end};
Q=Q{end};