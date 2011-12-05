function R=randmio_und_signed(R, ITER)
%RANDMIO_UND_SIGNED     Random graph with preserved degree distribution
%
%   R = randmio_und_signed(A,ITER);
%
%   This function randomizes an undirected weighted network with positive
%   and negative weights, while simultaneously preserving the degree 
%   distribution of positive and negative weights. The function does not 
%   preserve the strength distribution in weighted networks.
%
%   Input:      A,      undirected (binary/weighted) connection matrix
%               ITER,   rewiring parameter
%                       (each edge is rewired approximately ITER times)
%
%   Output:     R,      randomized network
%
%   Reference: Maslov and Sneppen (2002) Science 296:910
%
%
%   2011
%   Dani Bassett, UCSB
%   Mika Rubinov, UNSW

%   Modification History:
%   Mar 2011: Original (based on randmio_und.m)

[i j]=find(tril(R));
K=length(i);
[i_plus j_plus] = find(tril(R)>0);
[i_minus j_minus] = find(tril(R)<0);
K_plus = length(i_plus);
K_minus = length(i_minus);

ITER=K*ITER;

for iter=1:ITER
    while 1                                     %while not rewired
        while 1
            % choose two edges to rewire - but make sure they are either
            % both positive or both negative
            if rand>0.5 % chooses to rewire positive weighs and negative weights at equal rates
                e1=ceil(K_plus*rand);
                e2=ceil(K_plus*rand);
                type = 1;
            else
                e1=ceil(K_minus*rand);
                e2=ceil(K_minus*rand);
                type = 2;
            end
            if type==1;
                while (e2==e1),
                    e2=ceil(K_plus*rand);
                end
                a=i_plus(e1); b=j_plus(e1);
                c=i_plus(e2); d=j_plus(e2);
                if all(a~=[c d]) && all(b~=[c d]);
                    break           %all four vertices must be different
                end
            end
            if type==2;
                while (e2==e1),
                    e2=ceil(K_minus*rand);
                end
                a=i_minus(e1); b=j_minus(e1);
                c=i_minus(e2); d=j_minus(e2);        
                if all(a~=[c d]) && all(b~=[c d]);
                    break           %all four vertices must be different
                end
            end
        end
        if type==1;
            if rand>0.5
                i_plus(e2)=d; j_plus(e2)=c; 	%flip edge c-d with 50% probability
                c=i_plus(e2); d=j_plus(e2); 	%to explore all potential rewirings
            end
            %rewiring condition
            if ~(R(a,d) || R(c,b))
                R(a,d)=R(a,b); R(a,b)=0;
                R(d,a)=R(b,a); R(b,a)=0;
                R(c,b)=R(c,d); R(c,d)=0;
                R(b,c)=R(d,c); R(d,c)=0;
                
                j_plus(e1) = d;          %reassign edge indices
                j_plus(e2) = b;
                break;
            end %rewiring condition
        end
        if type==2;
            if rand>0.5
                i_minus(e2)=d; j_minus(e2)=c; 	%flip edge c-d with 50% probability
                c=i_minus(e2); d=j_minus(e2); 	%to explore all potential rewirings
            end
            %rewiring condition
            if ~(R(a,d) || R(c,b))
                R(a,d)=R(a,b); R(a,b)=0;
                R(d,a)=R(b,a); R(b,a)=0;
                R(c,b)=R(c,d); R(c,d)=0;
                R(b,c)=R(d,c); R(d,c)=0;
                
                j_minus(e1) = d;          %reassign edge indices
                j_minus(e2) = b;
                break;
            end %rewiring condition
        end
    end %while not rewired
end %iterations