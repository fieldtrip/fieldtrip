function R=randmio_und(R, ITER)
%RANDMIO_UND     Random graph with preserved degree distribution
%
%   R = randmio_und(A,ITER);
%
%   This function randomizes an undirected network, while preserving the 
%   degree distribution. The function does not preserve the strength 
%   distribution in weighted networks.
%
%   Input:      A,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   2007-2011
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Apr 2008: Edge c-d is flipped with 50% probability, allowing to explore
%             all potential rewirings (Jonathan Power)


[i j]=find(tril(R));
K=length(i);
ITER=K*ITER;

for iter=1:ITER
    while 1                                     %while not rewired
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
            break;
        end %rewiring condition
    end %while not rewired
end %iterations