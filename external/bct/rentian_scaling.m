function [N,E] = rentian_scaling(A,XYZ,n)
%RENTIAN_SCALING    Physical Rentian scaling
%
% [N E] = rentian_scaling(A,XYZ,n)
%
% Physical Rentian scaling (or more simply Rentian scaling) is a property 
% of systems that are cost-efficiently embedded into physical space. It is 
% what is called a "topo-physical" property because it combines information 
% regarding the topological organization of the graph with information 
% about the physical placement of connections. Rentian scaling is present 
% in very large scale integrated circuits, the C. elegans neuronal network, 
% and morphometric and diffusion-based graphs of human anatomical networks.
% Rentian scaling is determined by partitioning the system into cubes, 
% counting the number of nodes inside of each cube (N), and the number of 
% edges traversing the boundary of each cube (E). If the system displays 
% Rentian scaling, these two variables N and E will scale with one another 
% in loglog space. The Rent's exponent is given by the slope of log10(E) 
% vs. log10(N), and can be reported alone or can be compared to the 
% theoretical minimum Rent's exponent to determine how cost efficiently the 
% network has been embedded into physical space. Note: if a system displays 
% Rentian scaling, it does not automatically mean that the system is 
% cost-efficiently embedded (although it does suggest that). Validation 
% occurs when comparing to the theoretical minimum Rent's exponent for that
% system.
%
% Inputs:
% 	A       MxM adjacency matrix 
%           must be unweighted, binary, and symmetric.
% 	XYZ     Vector of node placement coordinates
%           must be Mx3 matrix, where M is the number of nodes.
% 	n       Number of partitions to compute. Each partition is a data
%           point. You want a large enough number to adequately estimate
%           the Rent's exponent.
%
% Outputs:
%	N       nx1 vector of the number of nodes in each of the n partitions.
%	E       nx1 vector of the number of edges crossing the boundary of each
%           partition.
%
% Subsequent Analysis:
%   Rentian scaling plots are then created by: figure; loglog(E,N,'*');
%
%	To determine the Rent's exponent, p, it is important not to use
%	partitions which may be affected by boundary conditions. In Bassett et
%	al. 2010 PLoS CB, only partitions with N<M/2 were used in the
%	estimation of the Rent's exponent. Thus, we can define N_prime =
%	N(find(N<M/2)) and E_prime = E(find(N<M/2)). Next we need to determine
%	the slope of Eprime vs. Nprime in loglog space, which is the Rent's
%	exponent. There are many ways of doing this with more or less
%	statistical rigor. Robustfit in MATLAB is one such option:
%       [b,stats] = robustfit(log10(N_prime),log10(E_prime))
%   
%   Then the Rent's exponent is b(1,2) and the standard error of the
%   estimation is given by stats.se(1,2). 
%
% Note: n=5000 was used in Bassett et al. 2010 in PLoS CB.
%
%
% Reference: 
%   Danielle S. Bassett, Daniel L. Greenfield, Andreas Meyer-Lindenberg,
%   Daniel R. Weinberger, Simon W. Moore, Edward T. Bullmore. Efficient
%   physical embedding of topologically complex information processing
%   networks in brains and computer circuits. PLoS Comput Biol, 2010,
%   6(4):e1000748.
%
%
% Danielle Bassett, UCSB, 2010


% determine the number of nodes in the system
M = numel(XYZ(:,1)); 
% rescale coordinates so that they are all greater than unity
XYZn = XYZ-repmat(min(XYZ)-1,M,1);
% find the absolute minimum and maximum over all directions
nmax = max(max(XYZn));
nmin = min(min(XYZn));

% initialize variables
count = 0;
N = zeros(n,1); 
E = zeros(n,1);

% create partitions, and count the number of nodes inside the partition (N) and the number of % edges traversing the boundary of the partition  (E)
while count<(n+1);
    % define cube end points
    randx = sort((1+nmax-nmin).*rand(2,1),'ascend');
    % find nodes in cube
    L = find(XYZn(:,1)>randx(1) & XYZn(:,1)<randx(2) & XYZn(:,2)>randx(1) & XYZn(:,2)<randx(2) & XYZn(:,3)>randx(1) & XYZn(:,3)<randx(2));
    if ~isempty(L)
        count = count+1;
        % count edges crossing the boundary of the cube
        E(count,1) = sum(sum(A(L,setdiff(1:M,L))));
        % count nodes inside of the cube
        N(count,1) = numel(L);
    end
end

   