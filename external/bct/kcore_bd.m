function [CIJkcore,kn,peelorder,peellevel] = kcore_bd(CIJ,k)
%KCORE_BD       K-core
%
%   [CIJkcore,kn,peelorder,peellevel] = kcore_bd(CIJ,k);
%
%   The k-core is the largest subnetwork comprising nodes of degree at
%   least k. This function computes the k-core for a given binary directed
%   connection matrix by recursively peeling off nodes with degree lower
%   than k, until no such nodes remain.
%
%   input:          CIJ,        connection/adjacency matrix (binary, directed)
%                     k,        level of k-core
%
%   output:    CIJkcore,        connection matrix of the k-core.  This matrix
%                               only contains nodes of degree at least k.
%                    kn,        size of k-core
%                    peelorder, indices in the order in which they were
%                               peeled away during k-core decomposition
%                    peellevel, corresponding level - nodes at the same
%                               level have been peeled away at the same time
%
%   'peelorder' and 'peellevel' are similar the the k-core sub-shells
%   described in Modha and Singh (2010).
%
%   References: e.g. Hagmann et al. (2008) PLoS Biology
%
%   Olaf Sporns, Indiana University, 2007/2008/2010/2012

%#ok<*ASGLU>
%#ok<*AGROW>

peelorder = [];
peellevel = [];
iter = 0;

while 1 

    % get degrees of matrix
    [id,od,deg] = degrees_dir(CIJ); 

    % find nodes with degree <k
    ff = find((deg<k)&(deg>0));
    
    % if none found -> stop
    if (isempty(ff)) break; end;            %#ok<SEPEX>

    % peel away found nodes
    iter = iter+1;
    CIJ(ff,:) = 0;
    CIJ(:,ff) = 0;
    
    peelorder = [peelorder; ff']; 
    peellevel = [peellevel; iter.*ones(1,length(ff))'];
    
end;

CIJkcore = CIJ;
kn = sum(deg>0);

