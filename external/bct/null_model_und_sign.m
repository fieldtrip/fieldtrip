function W0 = null_model_und_sign(W,ITER)
%NULL_MODEL_UND_SIGN     Random graphs with preserved weight, degree and
%                        strength distributions
%
%   W0 = null_model_und_sign(W,ITER);
%
%   This function randomizes an undirected network with positive and
%   negative weights, while preserving the degree and strength
%   distributions. This function calls randmio_und.m
%
%   Input:      W,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     W0,     randomized network
%
%   References: Rubinov and Sporns (2011) NeuroImage
%
%
%   2011, Mika Rubinov, UNSW

%   Modification History
%   Mar 2011: Original


n=size(W,1);                                        %number of nodes
W(1:n+1:end)=0;                                     %clear diagonal

Ap=W>0;                                             %positive adjacency matrix
Ap_r=randmio_und(Ap,ITER);                          %randomized Ap
An=~Ap; An(1:n+1:end)=0;                            %negative adjacency matrix
An_r=~Ap_r; An_r(1:n+1:end)=0;                      %randomized negative adjacency

W0=zeros(n);                                        %null model network
for s=[1 -1]                                        
    switch s                                        %switch sign (positive/negative)
        case 1
            S=sum(W.*Ap);                           %positive strength
            Wv=sort(W(triu(Ap)));                   %ordered weights vector
            [I J]=find(triu(Ap_r));                 %weights indices
            Lij=n*(J-1)+I;                          %linear weights indices
        case -1
            S=sum(-W.*An);                          %negative strength
            Wv=sort(-W(triu(An)));                  %ordered weights vector
            [I J]=find(triu(An_r));                 %weights indices
            Lij=n*(J-1)+I;                          %linear weights indices
    end

    P=(S.'*S);                                      %expected weights matrix
    for m=numel(Wv):-1:1                            %iteratively explore all weights
        [dum Oind]=sort(P(Lij));                    %get indices of Lij that sort P
        r=ceil(rand*m);                             
        o=Oind(r);                                  %choose random index of sorted expected weight
        W0(Lij(o)) = s*Wv(r);                       %assign corresponding sorted weight at this index

        f = 1 - Wv(r)/S(I(o));                      %readjust expected weight probabilities for node I(o)
        P(I(o),:) = P(I(o),:)*f;                    %[1 - Wv(r)/S(I(o)) = (S(I(o)) - Wv(r))/S(I(o))]
        P(:,I(o)) = P(:,I(o))*f;
        f = 1 - Wv(r)/S(J(o));                      %readjust expected weight probabilities for node J(o)
        P(J(o),:) = P(J(o),:)*f;                    %[1 - Wv(r)/S(J(o)) = (S(J(o)) - Wv(r))/S(J(o))]
        P(:,J(o)) = P(:,J(o))*f;                    
        
        S([I(o) J(o)]) = S([I(o) J(o)])-Wv(r);      %readjust strengths of nodes I(o) and J(o)
        Lij(o)=[];                                  %remove current index from further consideration
        I(o)=[];
        J(o)=[];
        Wv(r)=[];                                   %remove current weight from further consideration
    end
end
W0=W0+W0.';
