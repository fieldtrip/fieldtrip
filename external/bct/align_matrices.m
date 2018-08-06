function [Mreordered,Mindices,cost] = align_matrices(M1,M2,dfun,flag)
%ALIGN_MATRICES         alignment of two matrices
%
%   [Mreordered,Mindices,cost] = align_matrices(M1,M2,dfun,flag)
%
%   This function aligns two matrices relative to one another by reordering
%   the nodes in M2.  The function uses a version of simulated annealing.
%
%   Inputs:     M1             = first connection matrix (square)
%               M2             = second connection matrix (square)
%               dfun           = distance metric to use for matching:
%                                'absdff' = absolute difference
%                                'sqrdff' = squared difference
%                                'cosang' = cosine of vector angle
%
%               Mreordered     = reordered connection matrix M2
%               Mindices       = reordered indices
%               cost           = distance between M1 and Mreordered
%
%   Connection matrices can be weighted or binary, directed or undirected.
%   They must have the same number of nodes.  M1 can be entered in any
%   node ordering.
%
%   Note that in general, the outcome will depend on the initial condition
%   (the setting of the random number seed).  Also, there is no good way to 
%   determine optimal annealing parameters in advance - these parameters 
%   will need to be adjusted "by hand" (particularly H, Texp, T0, and Hbrk).  
%   For large and/or dense matrices, it is highly recommended to perform 
%   exploratory runs varying the settings of 'H' and 'Texp' and then select 
%   the best values.
%
%   Based on extensive testing, it appears that T0 and Hbrk can remain
%   unchanged in most cases.  Texp may be varied from 1-1/H to 1-10/H, for
%   example.  H is the most important parameter - set to larger values as
%   the problem size increases.  Good solutions can be obtained for
%   matrices up to about 100 nodes.  It is advisable to run this function
%   multiple times and select the solution(s) with the lowest 'cost'.
%
%   If the two matrices are related it may be very helpful to pre-align them
%   by reordering along their largest eigenvectors:
%       [v,~] = eig(M1); v1 = abs(v(:,end)); [a1,b1] = sort(v1);
%       [v,~] = eig(M2); v2 = abs(v(:,end)); [a2,b2] = sort(v2);
%       [a,b,c] = overlapMAT2(M1(b1,b1),M2(b2,b2),'dfun',1);
%
%   Setting 'Texp' to zero cancels annealing and uses a greedy algorithm
%   instead.
%
%   Yusuke Adachi, University of Tokyo, 2010
%   Olaf Sporns, Indiana University, 2010

N = size(M1,1);

% define maxcost (greatest possible difference)
switch dfun
case 'absdff'
    maxcost = sum(abs(sort(M1(:))-(sort(M2(:),'descend'))));
case 'sqrdff'
    maxcost = sum((sort(M1(:))-(sort(M2(:),'descend'))).^2);
case 'cosang'
    maxcost = pi/2;
end;

% initialize lowcost
switch dfun
case 'absdff'
    lowcost = sum(sum(abs(M1-M2)))/maxcost;
case 'sqrdff'
    lowcost = sum(sum((M1-M2).^2))/maxcost;
case 'cosang'
    lowcost = acos(dot(M1(:),M2(:))./sqrt(dot(M1(:),M1(:))*dot(M2(:),M2(:))))/maxcost;
end;

% initialize 
mincost = lowcost;
anew = 1:N;
amin = 1:N;
h = 0; hcnt = 0;

% set annealing parameters
% H determines the maximal number of steps
% Texp determines the steepness of the temperature gradient
% T0 sets the initial temperature (and scales the energy term)
% Hbrk sets a break point for the simulation (no further improvement)
H = 1e06; Texp = 1-1/H; T0 = 1e-03; Hbrk = H/10;
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
    switch dfun
        case 'absdff'
            costnew = sum(sum(abs(M1-M2(atmp,atmp))))/maxcost;
        case 'sqrdff'
            costnew = sum(sum((M1-M2(atmp,atmp)).^2))/maxcost;
        case 'cosang'
            M2atmp = M2(atmp,atmp);
            costnew = acos(dot(M1(:),M2atmp(:))./sqrt(dot(M1(:),M1(:))*dot(M2atmp(:),M2atmp(:))))/maxcost;
    end;
    % annealing step
    if (costnew < lowcost) || (rand < exp(-(costnew-lowcost)/T))
        anew = atmp;
        lowcost = costnew;
        % is this the absolute best?
        if (lowcost<mincost)
            amin = anew;
            mincost = lowcost;
            if (flag==1) 
                disp(['step ',num2str(h),' ... current lowest cost = ',num2str(mincost)]);
            end;
            hcnt = 0;
        end;
        % if the cost is 0 we're done
        if (mincost==0)
            break;
        end;
    end;
end;
disp(['step ',num2str(h),' ... final lowest cost = ',num2str(mincost)]);

% prepare output
Mreordered = M2(amin,amin);
Mindices = amin;
cost = mincost;

