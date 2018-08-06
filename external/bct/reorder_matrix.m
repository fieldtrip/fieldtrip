function [Mreordered,Mindices,cost] = reorder_matrix(M1,cost,flag)
%REORDER_MATRIX         matrix reordering for visualization
%
%   [Mreordered,Mindices,cost] = reorder_matrix(M1,cost,flag)
%
%   This function rearranges the nodes in matrix M1 such that the matrix
%   elements are squeezed along the main diagonal.  The function uses a
%   version of simulated annealing. 
%
%   Inputs:     M1             = connection matrix (weighted or binary, 
%                                directed or undirected)
%               cost           = 'line' or 'circ', for shape of lattice
%                                cost (linear or ring lattice)
%
%               Mreordered     = reordered connection matrix
%               Mindices       = reordered indices
%               cost           = distance between M1 and Mreordered
%
%   Note that in general, the outcome will depend on the initial condition
%   (the setting of the random number seed).  Also, there is no good way to 
%   determine optimal annealing parameters in advance - these paramters 
%   will need to be adjusted "by hand" (particularly H, Texp, and T0).  
%   For large and/or dense matrices, it is highly recommended to perform 
%   exploratory runs varying the settings of 'H' and 'Texp' and then select 
%   the best values.
%
%   Based on extensive testing, it appears that T0 and Hbrk can remain
%   unchanged in most cases.  Texp may be varied from 1-1/H to 1-10/H, for
%   example.  H is the most important parameter - set to larger values as
%   the problem size increases.  It is advisable to run this function
%   multiple times and select the solution(s) with the lowest 'cost'.
%
%   Setting 'Texp' to zero cancels annealing and uses a greedy algorithm
%   instead.
%
%   Yusuke Adachi, University of Tokyo 2010
%   Olaf Sporns, Indiana University 2010

N = size(M1,1);

% generate cost function
if (strcmp(cost,'line'))
    profil = fliplr(normpdf(1:N,0,N/2));
end;
if (strcmp(cost,'circ'))
    profil = fliplr(normpdf(1:N,N/2,N/4));
end;
COST = (toeplitz(profil,profil).*~eye(N));
COST = COST./sum(sum(COST));

% establish maxcost, lowcost, mincost
maxcost = sum(sort(COST(:)).*(sort(M1(:))));
lowcost = sum(sum(M1.*COST))/maxcost;
mincost = lowcost;

% initialize
anew = 1:N;
amin = 1:N;
h = 0; hcnt = 0;

% set annealing parameters
% H determines the maximal number of steps
% Texp determines the steepness of the temperature gradient
% T0 sets the initial temperature (and scales the energy term)
% Hbrk sets a break point for the simulation (if no further improvement)
H = 1e04; Texp = 1-10/H; T0 = 1e-03; Hbrk = H/10;
%Texp = 0;

while h<H
    h = h+1; hcnt = hcnt+1;
    % terminate if no new mincost has been found for some time
    if (hcnt>Hbrk)
        break; 
    end;
    % current temperature
    T = T0*Texp^h;
    % choose two positions at random and flip them
    atmp = anew;
    %r = randperm(N);  % slower
    r = ceil(rand(1,2).*N);
    atmp(r(1)) = anew(r(2));
    atmp(r(2)) = anew(r(1));
    costnew = sum(sum(M1(atmp,atmp).*COST))/maxcost;
    % annealing
    if (costnew < lowcost) || (rand < exp(-(costnew-lowcost)/T))
        anew = atmp;
        lowcost = costnew;
        % is this a new absolute best?
        if (lowcost<mincost)
            amin = anew;
            mincost = lowcost;
            if (flag==1) 
                disp(['step ',num2str(h),' ... current lowest cost = ',num2str(mincost)]);
            end;
            hcnt = 0;
        end;
    end;
end;
disp(['step ',num2str(h),' ... final lowest cost = ',num2str(mincost)]);

% prepare output
Mreordered = M1(amin,amin);
Mindices = amin;
cost = mincost;

