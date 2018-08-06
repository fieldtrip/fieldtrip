function [B,E,K] = evaluate_generative_model(A,Atgt,D,modeltype,modelvar,params)
% EVALUATE_GENERATIVE_MODEL     generate and evaluate synthetic networks
%
%   [B,E,K] = EVALUATE_GENERATIVE_MODEL(A,Atgt,D,m,modeltype,modelvar,params) 
%
%   Generates synthetic networks and evaluates their energy function (see
%   below) using the models described in the study by Betzel et al (2016)
%   in Neuroimage.
%
%   Inputs:
%           A,          binary network of seed connections
%           Atgt,       binary network against which synthetic networks are
%                       compared
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
%
%   Outputs:
%           B,          m x number of networks matrix of connections
%           E,          energy for each synthetic network
%           K,          Kolmogorov-Smirnov statistics for each synthetic
%                       network.
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
%   Note: Energy is calculated in exactly the same way as in Betzel et
%   al (2016). There are four components to the energy are KS statistics
%   comparing degree, clustering coefficient, betweenness centrality, and 
%   edge length distributions. Energy is calculated as the maximum across
%   all four statistics.
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015

m = nnz(Atgt)/2;
n = length(Atgt);
x = cell(4,1);
x{1} = sum(Atgt,2);
x{2} = clustering_coef_bu(Atgt);
x{3} = betweenness_bin(Atgt)';
x{4} = D(triu(Atgt,1) > 0);

B = generative_model(A,D,m,modeltype,modelvar,params);
nB = size(B,2);

K = zeros(nB,4);
for iB = 1:nB
    b = zeros(n);
    b(B(:,iB)) = 1;
    b = b + b';
    y = cell(4,1);
    y{1} = sum(b,2);
    y{2} = clustering_coef_bu(b);
    y{3} = betweenness_bin(b)';
    y{4} = D(triu(b,1) > 0);
    for j = 1:4
        K(iB,j) = fcn_ks(x{j},y{j});
    end
end
E = max(K,[],2);


function kstat = fcn_ks(x1,x2)
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);

