function Z=module_degree_zscore(W,Ci,flag)
%MODULE_DEGREE_ZSCORE       Within-module degree z-score
%
%   Z=module_degree_zscore(W,Ci,flag);
%
%   The within-module degree z-score is a within-module version of degree
%   centrality.
%
%   Inputs:     W,      binary/weighted, directed/undirected connection matrix
%               Ci,     community affiliation vector
%               flag,   0, undirected graph (default)
%                       1, directed graph: out-degree
%                       2, directed graph: in-degree
%                       3, directed graph: out-degree and in-degree
%
%   Output:     Z,      within-module degree z-score.
%
%   Reference: Guimera R, Amaral L. Nature (2005) 433:895-900.
%
%
%   Mika Rubinov, UNSW, 2008-2010

if ~exist('flag','var')
    flag=0;
end

switch flag
    case 0; % no action required
    case 1; % no action required
    case 2; W=W.';
    case 3; W=W+W.';
end

n=length(W);                        %number of vertices
Z=zeros(n,1);
for i=1:max(Ci)
    Koi=sum(W(Ci==i,Ci==i),2);
    Z(Ci==i)=(Koi-mean(Koi))./std(Koi);
end

Z(isnan(Z))=0;