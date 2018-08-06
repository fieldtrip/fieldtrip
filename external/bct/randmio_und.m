function [R,eff]=randmio_und(R, ITER)
%RANDMIO_UND     Random graph with preserved degree distribution
%
%   R = randmio_und(W,ITER);
%   [R eff]=randmio_und(W, ITER);
%
%   This function randomizes an undirected network, while preserving the 
%   degree distribution. The function does not preserve the strength 
%   distribution in weighted networks.
%
%   Input:      W,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%               eff,    number of actual rewirings carried out
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   2007-2012
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Olaf Sporns, IU

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Apr 2008: Edge c-d is flipped with 50% probability, allowing to explore
%             all potential rewirings (Jonathan Power)
%   Mar 2012: Limit number of rewiring attempts, count number of successful
%             rewirings (Olaf Sporns)


n=size(R,1);
[i,j]=find(tril(R));
K=length(i);
ITER=K*ITER;

% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)));
% actual number of successful rewirings
eff = 0;

for iter=1:ITER
    att=0;
    while (att<=maxAttempts)                                     %while not rewired
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1),
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);

            if all(a~=[c d]) && all(b~=[c d]);
                break           %all four vertices must be different
            end
        end

        if rand>0.5
            i(e2)=d; j(e2)=c; 	%flip edge c-d with 50% probability
            c=i(e2); d=j(e2); 	%to explore all potential rewirings
        end
        
        %rewiring condition
        if ~(R(a,d) || R(c,b))
            R(a,d)=R(a,b); R(a,b)=0;
            R(d,a)=R(b,a); R(b,a)=0;
            R(c,b)=R(c,d); R(c,d)=0;
            R(b,c)=R(d,c); R(d,c)=0;

            j(e1) = d;          %reassign edge indices
            j(e2) = b;
            eff = eff+1;
            break;
        end %rewiring condition
        att=att+1;
    end %while not rewired
end %iterations